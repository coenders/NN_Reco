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

from icecube import dataclasses, dataio, icetray
from icecube.phys_services import I3Calculator
from icecube.icetray import I3Units
import icecube.MuonGun
import numpy as np


def calc_depositedE(physics_frame):
    I3Tree = physics_frame['I3MCTree']
    losses = 0
    for p in I3Tree:
        if not p.is_cascade: continue
        if not p.location_type == dataclasses.I3Particle.InIce: continue
        if p.shape == p.Dark: continue
        if p.type in [p.Hadrons, p.PiPlus, p.PiMinus, p.NuclInt]:
            if p.energy < 1*I3Units.GeV:
                losses += 0.8*p.energy
            else:
                energyScalingFactor = 1.0 + ((p.energy/I3Units.GeV/0.399)**-0.130)*(0.467 - 1)
                losses += energyScalingFactor*p.energy
        else:
            losses += p.energy
    return losses


# TODO: Check if particle was in the detector volume at all
def classificationTag_MCTree(frame):
    def within_icecube(i3p):
        distance = I3Calculator.distance_along_track(i3p, origin)

        return distance - dist2icecube < i3p.length

    dist2icecube = 1e3
    min_track_length = 5
    origin = dataclasses.I3Position(0, 0, 0)

    # return codes
    CASCADE = 1
    TRACK = 2
    DOUBLE_BANG = 3

    # Analyze MCTree
    mctree = frame["I3MCTree"]

    # Get all neutrinos
    neutrinos = [p_i for p_i in mctree.get_primaries()
            if abs(p_i.pdg_encoding) in [12, 14, 16]]

    if len(neutrinos) < 1:
        raise IndexError("No Neutrinos in MCTree")

    # get neutrino with highest energy, should there be more than one, though?
    nu = neutrinos[0]
    for nu_i in neutrinos[1:]:
        if nu_i.energy > nu.energy:
            nu = nu_i

    # choose first particle close to IceCube (avoid regeneration of nutau)
    while abs(nu.pdg_encoding) in (12, 14, 16) or not within_icecube(nu):
        daughters = [d_i for d_i in mctree.get_daughters(nu)
                if d_i.shape != dataclasses.I3Particle.Dark]

        if len(daughters) < 1:
            return CASCADE

        nu = daughters[0]

    print "Neutrino", nu

    # search for particle with minimum feasible length
    while nu.length < min_track_length:
        daughters = [d_i for d_i in mctree.get_daughters(nu)
                if d_i.shape != dataclasses.I3Particle.Dark]

        if len(daughters) < 1:
            return CASCADE

        nu = daughters[0]

    print "Shapy", nu

    # check whether the particle exits IceCube
    endpoint = nu.pos + nu.length * nu.dir
    r = np.sqrt(endpoint.x**2 + endpoint.y**2)
    z = abs(endpoint.z)

    if max(z, r) > dist2icecube:
        # Track is leaving the building
        return TRACK

    # checking what happens with the "thing" inside of IceCube
    if abs(nu.pdg_encoding) == 13:
        # muon, nothing special expected to happen
        return TRACK
    elif abs(nu.pdg_encoding) == 15:
        # tau, check the decay process
        daughters = [d_i for d_i in mctree.get_daughters(nu)
                if d_i.shape != dataclasses.I3Particle.Dark]
        for d in daughters:
            if d.length > min_track_length:
                return Track
        return DOUBLE_BANG
    else:
        raise ValueError("Don't know to handle this particle of type",
                nu.pdg_encoding)


def classificationTag(physics_frame):
    energy = calc_depositedE(physics_frame)
    ParticelList = [12, 14, 16]
    I3Tree = physics_frame['I3MCTree']
    primary_list = I3Tree.get_primaries()
    if len(primary_list) == 1:
        neutrino = I3Tree[0]
    else:
        for p in primary_list:
            pdg = p.pdg_encoding
            if abs(pdg) in ParticelList:
                neutrino = p
    if abs(neutrino.pdg_encoding) == 12: # primary particle is a electron neutrino
        classificationTag = 1
    elif abs(neutrino.pdg_encoding) == 14: # primary particle is a muon neutrino
        classificationTag = 2
    elif abs(neutrino.pdg_encoding) == 16: # primary particle is a tau neutrino
        listChildren = I3Tree.children(I3Tree.first_child(neutrino))
        if not listChildren: # without this, the function collapses
            if energy > 10**6: # more than 1 PeV, due to energy in GeV
                classificationTag = 3
            else:
                classificationTag = 1
        else:
            for i in listChildren:
                if abs(i.pdg_encoding) == 13:
                     classificationTag = 2
                else:
                    if energy > 10**6: # more than 1 PeV, due to energy in GeV
                        classificationTag = 3
                    else:
                        classificationTag = 1
    else:
        print "Error: primary particle wasnt a neutrino"
    # classificationTag = 1 means cascade
    # classificationTag = 2 means track
    # classificationTag = 3 means double bang
    return classificationTag


def starting(physics_frame):
    gcdfile = "/data/sim/sim-new/downloads/GCD/GeoCalibDetectorStatus_2012.56063_V0.i3.gz"
    N = 0
    ParticelList = [12, 14, 16]
    I3Tree = physics_frame['I3MCTree']
    primary_list = I3Tree.get_primaries()
    if len(primary_list) == 1:
        neutrino = I3Tree[0]
    else:
        for p in primary_list:
            pdg = p.pdg_encoding
            if abs(pdg) in ParticelList:
                neutrino = p
    surface = icecube.MuonGun.ExtrudedPolygon.from_file(gcdfile, padding=-N)
    intersections = surface.intersection(neutrino.pos + neutrino.length*neutrino.dir, neutrino.dir)
    if intersections.first <= 0 and intersections.second > 0:
        starting = 0 # starting event
    else:
        starting = 1 # through-going or stopping event
    return starting

def up_or_down(physics_frame):
    zenith = physics_frame["LineFit"].dir.zenith
    if zenith > 1.5*np.pi or zenith < 0.5*np.pi:
        up_or_down = 1 # down-going
    else:
        up_or_down = 0 # up-going
    return up_or_down


def coincidenceLabel(physics_frame):
    primary_list = physics_frame["I3MCTree"].get_primaries()
    if len(primary_list) > 1:
        coincidence = 1
    else:
        coincidence = 0
    return coincidence




