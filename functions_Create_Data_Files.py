#!/usr/bin/env python
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
# may some packages are not needed


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
