# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 11:12:10 2023

@author: alexandre
"""

import pickle     as pk
import cblind     as cb
import MDAnalysis as mda
import numpy      as np
import matplotlib.pyplot as plt

from optparse import OptionParser

class Dens:
    
    def __init__(self, BINS, DENSITY):
        
        self.BINS    = BINS
        self.DENSITY = DENSITY
        

parser = OptionParser(usage       = 'python insertion.py -t <input .tpr/.gro/.pdb file> -x <input .xtc/.trr file> -o <output .png/.jpg file',
                      prog        = 'Insertion',
                      description = 'This program is used calculate the lateral density.')

parser.add_option('-t', '--tpr',
                  action = 'store', type = 'string',
                  default = '5us.tpr',
                  help = 'Input .tpr File')

parser.add_option('-x', '--xtc',
                  action = 'store', type = 'string',
                  default = '5us.xtc',
                  help = 'Input .xtc File')

parser.add_option('-o', '--outfile',
                  action = 'store', type = 'string',
                  default = 'density.png',
                  help = 'Output .png file')

parser.add_option('-i', '--initial',
                  action = 'store', type = int,
                  default = 0,
                  help = 'Initial time in ps')

parser.add_option('-f', '--final',
                  action = 'store', type = int,
                  default = 5000000,
                  help = 'Final time in ps')


options, args = parser.parse_args()


def Density(POS, MASS, BINS, V):

    OUT = []

    for i in range(BINS.shape[0]-1):

        idx_bot = POS >= BINS[i]
        idx_top = POS <  BINS[i+1]

        idx = idx_bot*idx_top

        masses = np.sum(MASS[idx])

        OUT.append(masses/V)

    return np.array(OUT)



if __name__ == '__main__':


    U = mda.Universe(options.tpr, options.xtc,
                     in_memory = False)

    PMB = U.select_atoms('resname PMB1',
                         periodic = True, updating = True)

    MEMBRANE = U.select_atoms('resname *EMP or resname POP* or resname CLD2',
                              periodic = True, updating = True)

    REMP_H  = U.select_atoms('resname REMP and (name GM* or name PO* or name S*)',
                             periodic = True, updating = True)

    REMP_T  = U.select_atoms('resname REMP and (name GL* or name C*)',
                             periodic = True, updating = True)
    
    FEMP_H  = U.select_atoms('resname FEMP and (name GM* or name PO* or name S*)',
                             periodic = True, updating = True)

    FEMP_T  = U.select_atoms('resname FEMP and (name GL* or name C*)',
                             periodic = True, updating = True)

    
    LIPID   = U.select_atoms('resname POP* or resname CDL2',
                             periodic = True, updating = True)

    idx_i =int(np.ceil(options.initial/U.trajectory.dt))
    idx_f =int(np.ceil(options.final/U.trajectory.dt))

    DIMS = []
    for frame in range(idx_i, idx_f, 1):
        U.trajectory[frame]

        MEMBRANE.universe.trajectory[frame]

        COM = MEMBRANE.center_of_mass()
        dims = U.dimensions[:3].copy() - COM
        DIMS.append([dims, -1*COM])

    DIMS = np.array(DIMS)
    MinD, MaxD = np.min(DIMS[:,1], axis = 0), np.max(DIMS[:,0], axis = 0)

    BINS = np.arange(MinD[2], MaxD[2] + 1, 1)

    DENS = []

    n_frames = idx_f - idx_i

    for f, frame in enumerate(U.trajectory[idx_i:idx_f]):

        print(f'Done {int(100*f/(n_frames-1))} %', end = '\r')

        x   = DIMS[f,0,0]
        y   = DIMS[f,0,1]
        z   = BINS[1]-BINS[0]
        V   = x*y*z
        com = DIMS[f,1,2]


        POS_PMB  = PMB.positions
        MASS_PMB = PMB.masses

        POS_REMPH  = REMP_H.positions
        MASS_REMPH = REMP_H.masses

        POS_REMPT  = REMP_T.positions
        MASS_REMPT = REMP_T.masses
        
        POS_FEMPH  = FEMP_H.positions
        MASS_FEMPH = FEMP_H.masses

        POS_FEMPT  = FEMP_T.positions
        MASS_FEMPT = FEMP_T.masses

        POS_LIPID  = LIPID.positions
        MASS_LIPID = LIPID.masses

        dens_pmb   = Density(POS_PMB[:,2]   + com, MASS_PMB,   BINS, V)
        dens_remph = Density(POS_REMPH[:,2] + com, MASS_REMPH, BINS, V)
        dens_rempt = Density(POS_REMPT[:,2] + com, MASS_REMPT, BINS, V)
        dens_femph = Density(POS_FEMPH[:,2] + com, MASS_FEMPH, BINS, V)
        dens_fempt = Density(POS_FEMPT[:,2] + com, MASS_FEMPT, BINS, V)
        dens_lipid = Density(POS_LIPID[:,2] + com, MASS_LIPID, BINS, V)


        DENS.append([dens_pmb,
                     dens_remph,
                     dens_rempt,
                     dens_femph,
                     dens_fempt,
                     dens_lipid])

    DENS = Dens(0.5*(BINS[1:] + BINS[:-1]), np.array(DENS))
    
    with open('density.pk', 'wb') as f:
        pk.dump(DENS, f, protocol = 4)
    f.close()