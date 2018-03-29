#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT /data/user/tglauch/Software/combo/build
# coding: utf-8

"""This file is part of DeepIceLearning
DeepIceLearning is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


from icecube import dataio, icetray
from scipy.stats import moment, skew, kurtosis
import numpy as np
import math
import tables
import argparse
import os, sys
from configparser import ConfigParser
from reco_quantities import *
import cPickle as pickle

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset_config",
        help="main config file, user-specific",
        type=str, default='default.cfg')
    parser.add_argument(
        "--files",
        help="files to be processed",
        type=str, nargs="+", required=False)
    parser.add_argument(
        "--filelist",
        help="Path to a filelist to be processed",
        type=str, required=False)
    parser.add_argument(
        "--version",
        action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()

    return args


args = parseArguments().__dict__

dataset_configparser = ConfigParser()
try:
    dataset_configparser.read(args['dataset_config'])
except Exception:
    raise Exception('Config File is missing!!!!')

# File paths #########
basepath = str(dataset_configparser.get('Basics', 'MC_path'))
geometry_file = str(dataset_configparser.get('Basics', 'geometry_file'))
outfolder = str(dataset_configparser.get('Basics', 'out_folder'))
pulsemap_key = str(dataset_configparser.get('Basics', 'PulseSeriesMap'))


def read_variables(cfg_parser):
    """Function reading a config file, defining the variables to be read
       from the MC files.

    Arguments:
    cfg_parser: config parser object for the config file

    Returns:
    dtype : the dtype object defining the shape and names of the MC output
    data_source: list defining the types,names and ranges of monte carlo data
                to be saved from a physics frame
                (e.g [('variable',['MCMostEnergeticTrack'].energy, [1e2,1e9])])
    """
    dtype = []
    data_source = []
    for i, key in enumerate(cfg_parser.keys()):
        if key == 'DEFAULT' or key == 'Basics':
            continue
        cut = [-np.inf, np.inf]
        if 'min' in cfg_parser[key].keys():
            cut[0] = float(cfg_parser[key]['min'])
        if 'max' in cfg_parser[key].keys():
            cut[1] = float(cfg_parser[key]['max'])
        if 'variable' in cfg_parser[key].keys():
            data_source.append(('variable', cfg_parser[key]['variable'], cut))
        elif 'function' in cfg_parser[key].keys():
            data_source.append(('function', cfg_parser[key]['function'], cut))
        else:
            raise Exception(
                'No Input Type given. Variable or funtion must be given')
        dtype.append((str(key), eval('np.' + cfg_parser[key]['out_type'])))
    dtype = np.dtype(dtype)

    return dtype, data_source


def preprocess_grid(geometry):
    # rotate IC into x-y-plane
    dom_6_pos = geometry[icetray.OMKey(6, 1)].position
    dom_1_pos = geometry[icetray.OMKey(1, 1)].position
    theta = -np.arctan(
        (dom_6_pos.y - dom_1_pos.y) / (dom_6_pos.x - dom_1_pos.x))
    c, s = np.cos(theta), np.sin(theta)
    rot_mat = np.matrix([[c, -s], [s, c]])

    # om > 60 are icetops  om 79-87 are deepcore --> exclude
    DOM_List = sorted(
        [i for i in geometry.keys()
         if i.om < 61 and i.string not in range(79, 87)])
    xpos = [geometry[i].position.x for i in DOM_List]
    ypos = [geometry[i].position.y for i in DOM_List]
    zpos = [geometry[i].position.z for i in DOM_List]

    rotxy = [np.squeeze(np.asarray(np.dot(rot_mat, xy)))
             for xy in zip(xpos, ypos)]
    xpos, ypos = zip(*rotxy)
    return xpos, ypos, zpos, DOM_List


def make_grid_dict(input_shape, geometry):
    """Put the Icecube Geometry in a cuboid grid.
    For each DOM calculate the corresponding grid position.
    Rotates the x-y-plane in order to make icecube better fit into a grid.

    Arguments:
    input_shape : The shape of the grid (x,y,z)
    geometry : Geometry file containing the positions of the DOMs in
    the Detector

    Returns:
    grid: a dictionary mapping (string, om) => (grid_x, grid_y, grid_z),
    i.e. dom id to its index position in the cuboid grid
    dom_list_ret: list of all (string, om), i.e. list of dom ids in the geofile
    (dom_list_ret==sorted(grid.keys()))
    """
    grid = dict()
    xpos, ypos, zpos, DOM_List = preprocess_grid(geometry)

    xmin, xmax = np.min(xpos), np.max(xpos)
    delta_x = (xmax - xmin) / (input_shape[0] - 1)
    xmin, xmax = xmin - delta_x / 2, xmax + delta_x / 2
    ymin, ymax = np.min(ypos), np.max(ypos)
    delta_y = (ymax - ymin) / (input_shape[1] - 1)
    ymin, ymax = ymin - delta_y / 2, ymax + delta_y / 2
    zpos_reshaped = np.array(zpos).reshape(78, 60)
    zmin, zmax = np.median(map(np.min, zpos_reshaped)),
    np.median(map(np.max, zpos_reshaped))
    delta_z = (zmax - zmin) / (input_shape[2] - 1)
    zmin, zmax = zmin - delta_z / 2, zmax + delta_z / 2
    dom_list_ret = []
    for i, odom in enumerate(DOM_List):
        dom_list_ret.append((odom.string, odom.om))
# for all x,y,z-positions the according grid position is calculated and
# stored. the doms that lie outside the z-range are put in to the closest bin
# (see: https://www.dropbox.com/s/fsjuxrua28dz2fi/zbinning.png)
# z coordinates count from bottom to top (righthanded coordinate system)
        grid[(odom.string, odom.om)] = \
            (int(math.floor((xpos[i] - xmin) / delta_x)),
             int(math.floor((ypos[i] - ymin) / delta_y)),
             input_shape[2] - 1 - max(
                min(int(math.floor((zpos[i] - zmin) / delta_z)),
                    input_shape[2] - 1), 0))
    return grid, dom_list_ret


def make_autoHexGrid(geometry):
    """Put the Icecube Geometry in a rectangular grid.
    For each DOM calculate corresponding grid position. Rotates the x-y-plane
    in order to make icecube better fit into a grid.
    Method: aligns IC-strings which are not on the hexagonal grid + shifts
    x_positions such that no unfilled holes appear in the grid but rather empty
    edges (reduces dimensionality of the input and makes pattern recognition
    much easier)

    Arguments:
    geometry : Geometry file containing the
    positions of the DOMs in the Detector

    Returns:
    grid: a dictionary mapping (string, om) =>(grid_x, grid_y, grid_z),
    i.e. dom id to its index position in the cubic grid
    dom_list_ret: list of all (string, om),
    i.e. list of all dom ids in the geofile
    (sorted(dom_list_ret)==sorted(grid.keys()))
    """

    grid = dict()
    # assumes the standard IC shape:
    max_string = max(o.string for o in geometry.keys())
    max_dom = max(o.om for o in geometry.keys())
    if max_string < 78 or max_dom < 60:
        print "Define your own input_shape, makeHexGrid is only for standardIC"
        raise NameError('Wrong geometry file for standard IC processing')

    xpos, ypos, zpos, DOM_List = preprocess_grid(geometry)
    deltax = abs(xpos[0] - xpos[60])  # inserted by hand, any better idea ?
    deltay = abs(ypos[360] - ypos[0])

    nxRows, nyRows = 20, 10  # again, standard IC geometry (20x10 w/ holes)
    # align strings which do not lie on the hexagonal grid:
    xBands = np.linspace(np.amin(xpos) - deltax / 4.,
                         np.amax(xpos) + deltax / 4., nxRows + 1)
    yBands = np.linspace(np.amin(ypos) - deltay / 2.,
                         np.amax(ypos) + deltay / 2., nyRows + 1)
    xIndcs = np.digitize(xpos, xBands)
    yIndcs = np.digitize(ypos, yBands)
    # reset positions to the exact hex-grid positions
    xpos_aligned = deltax / 4. * xIndcs
    ypos_aligned = deltay / 2. * yIndcs

    # update deltas
    deltax_aligned = abs(xpos_aligned[0] - xpos_aligned[60])
    deltay_aligned = abs(ypos_aligned[360] - ypos_aligned[0])

    # shift the x-positions of each DOM to shift the hex-grid to a rect-grid
    xpos_shifted = xpos_aligned + deltax_aligned / 2. *\
        np.floor((ypos_aligned - (
            np.amin(ypos_aligned) + 1e-5)) / deltay_aligned)
    # center the new grid
    x_final = xpos_shifted - np.mean(xpos_shifted)
    y_final = ypos_aligned - np.mean(xpos_aligned)

    # final grid:
    xinput_bins = np.linspace(np.amin(x_final) - deltax_aligned / 2.,
                              np.amax(x_final) + deltax_aligned / 2.,
                              12)
    yinput_bins = np.linspace(np.amin(y_final) - deltay_aligned / 2.,
                              np.amax(y_final) + deltay_aligned / 2.,
                              11)
    zinput_bins = np.linspace(np.amin(zpos), np.amax(zpos), 60)

    dom_list_ret = []
    for i, odom in enumerate(DOM_List):
        dom_list_ret.append((odom.string, odom.om))
        grid[(odom.string, odom.om)] = \
            (np.digitize([x_final[i]], xinput_bins)[0],
             np.digitize([y_final[i]], yinput_bins)[0],
             np.digitize([zpos[i]], zinput_bins)[0])
    return grid, dom_list_ret


def analyze_grid(grid):
    """
    if you want to see which string/om the bins contain
    """
    dims = []
    for dim in range(3):
        for index in range(input_shape[dim]):
            strings = set()
            dims.append(list())
            for k, v in grid.items():
                if v[dim] == index:
                    if dim == 2:
                        strings.add(k[1])  # print om
                    else:
                        strings.add(k[0])  # print string
            dims[dim].append(strings)
    for i, c in enumerate("xyz"):
        print c
        for index, strings in enumerate(dims[i]):
            print index, strings


if __name__ == "__main__":

    # Raw print arguments
    print"\n ############################################"
    print("You are running the script with arguments: ")
    for a in args.keys():
        print(str(a) + ": " + str(args[a]))
    print"############################################\n "

    geo = dataio.I3File(geometry_file).pop_frame()['I3Geometry'].omgeo

    input_shape_par = dataset_configparser.get('Basics', 'input_shape')

    input_shape = [12, 11, 61]
    grid, DOM_List = make_autoHexGrid(geo)

    # Create HDF5 File ##########

    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    if args['filelist'] is not None :
        filelist = pickle.load(open(args['filelist'], 'r'))
        outfile = args['filelist'].replace('.pickle', '.h5')
    elif args['files'] is not None:
        filelist = args['files']
        tmp = filelist[0].replace('.i3.bz2', '.h5')
        outfile = outfolder+"/"+tmp.split("/")[-1] 
    else:
        raise Exception('No input files given')

    if os.path.exists(outfile):
        os.remove(outfile)

    dtype, data_source = read_variables(dataset_configparser)
    dtype_len = len(dtype)
    FILTERS = tables.Filters(complib='zlib', complevel=9)
    N_SLICES = 24
    T_MAX = 6000 # this cuts out some late pulses. But do we want/need them?
    time_bins = np.linspace(0, T_MAX, N_SLICES+1) ##try around with logspace?
                                                #More resolution for early pulses!
    print("Time-bins used for dicrete time-steps", time_bins)
    with tables.open_file(
        outfile, mode="w", title="Events for training the NN",
            filters=FILTERS) as h5file:
        charge = h5file.create_earray(
            h5file.root, 'charge', tables.Float64Atom(),
            (0, N_SLICES, input_shape[0], input_shape[1], input_shape[2], 1),
            title="Sum of charges per Dom per time-interval")
        reco_vals = tables.Table(h5file.root, 'reco_vals',
                                 description=dtype)
        h5file.root._v_attrs.shape = input_shape
        print('Created a new HDF File with the Settings:')
        print(h5file)

        np.save('grid.npy', grid)
        j = 0
        skipped_frames = 0
        for counter, f_name in enumerate(filelist):
            if counter % 10 == 0:
                print('Processing File {}/{}'.format(counter, len(filelist)))
            event_file = dataio.I3File(str(f_name), "r")
            print "Opening succesful"
            while event_file.more():
                physics_event = event_file.pop_physics()
                reco_arr = []
                for k, cur_var in enumerate(data_source):
                    if cur_var[0] == 'variable':
                        try:
                            cur_value = eval(
                                'physics_event{}'.format(cur_var[1]))
                        except Exception:
                            skipped_frames += 1
                            print('Attribute Error occured :{}'.
                                  format(cur_var[1]))
                            break

                    if cur_var[0] == 'function':
                        try:
                            cur_value = eval(
                                cur_var[1].replace('(x)', '(physics_event)'))
                        except Exception:
                            skipped_frames += 1
                            print(
                                'The given function is not implemented')
                            break

                    if cur_value < cur_var[2][0] or cur_value > cur_var[2][1]:
                        break
                    else:
                        reco_arr.append(cur_value)

                if not len(reco_arr) == dtype_len:
                    continue
                charge_arr = np.zeros(
                    (1, N_SLICES, input_shape[0], input_shape[1], input_shape[2], 1))
                 ###############################################
                pulses = physics_event[pulsemap_key].apply(physics_event)
                #times = []
                #for omkey in pulses.keys():
                #    for p in pulses[omkey]:
                #        times.append(p.time)
                #times_20pc = np.percentile(times,20)
                times_min = physics_event["SPEFit2TimeSplit1"].time #np.amin(times)

                #shift all times to the 20pc quantile and truncate at +1000ns
                #This could help to isolate the actual event
                for omkey in pulses.keys():
                    if (omkey.string, omkey.om) not in DOM_List:
                        continue
                    charges = np.zeros(N_SLICES+1)
                    for p in pulses[omkey]:
                        t_bin = np.digitize([p.time-times_min], time_bins)-1
                        charges[t_bin] += p.charge

                    gpos = grid[(omkey.string, omkey.om)]
                    for t_stamp, ch in enumerate(charges[:-1]):
                        charge_arr[0][t_stamp, gpos[0]][gpos[1]][gpos[2]] = ch

                charge.append(np.array(charge_arr))
                reco_vals.append(np.array(reco_arr))
                j += 1

            charge.flush()
            reco_vals.flush()
        print('###### Run Summary ###########')
        print('Processed: {} Frames \n Skipped {} \
            Frames with Attribute Error'.format(j, skipped_frames))
        h5file.root._v_attrs.len = j
        h5file.close()
