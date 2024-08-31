# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 17:03:48 2023

@author: alexandre
"""

import os
import numpy      as np
import pickle     as pk
import MDAnalysis as mda

from optparse import OptionParser

parser = OptionParser()
parser.add_option('-n', '--ndx_file',
                  default = 'POPC/Charmm/NDX/index_dns.ndx',
                  type    = str,
                  help    = 'input ndx file',
                  action  = 'store')
parser.add_option('-t', '--tpr_file',
                  default = 'POPC/Charmm/prod.tpr',
                  type    = str,
                  help    = 'input tpr file',
                  action  = 'store')
parser.add_option('-x', '--xtc_file',
                  default = 'POPC/Charmm/XTC/400_500.xtc',
                  type    = str,
                  help    = 'input xtc file',
                  action  = 'store')
parser.add_option('-s', '--selection',
                  default = 'Glycerol',
                  type    = str,
                  help    = 'selection group',
                  action  = 'store')
parser.add_option('-c', '--center',
                  default = 'PO4',
                  type    = str,
                  help    = 'centering group',
                  action  = 'store')
parser.add_option('-i', '--initial_time',
                  default = 400000,
                  type    = int,
                  help    = 'initial time',
                  action  = 'store')
parser.add_option('-f', '--final_time',
                  default = 410000,
                  type    = int,
                  help    = 'final time',
                  action  = 'store')
parser.add_option('-p', '--path',
                  default = './',
                  type    = str,
                  help    = 'master path',
                  action  = 'store')
parser.add_option('-o', '--outfile',
                  default = 'density.pk',
                  type    = str,
                  help    = 'file name',
                  action  = 'store')
options, arguments = parser.parse_args()

class Density:

    def __init__(self, X, DENSITY):

        self.X       = X
        self.DENSITY = DENSITY

def ReadNDX(path_to_file : str) -> dict:

    out = {}

    with open(path_to_file, 'r') as f:

        lines = f.readlines()

        line_number = 0

        while line_number < len(lines):

            line   = lines[line_number]
            fields = line.split()

            if fields[0] == '[':

                group      = fields[1]
                out[group] = []

                for l, line in enumerate(lines[line_number + 1:]):

                    if '[' in line:
                        line_number += (l + 1)
                        break

                    if line_number + (l+2) >= len(lines):
                        line_number = len(lines)

                    ids = line.split()
                    for _ in ids:
                        out[group].append(int(_) - 1)

    return out


def CalcBoxSize(U : mda.Universe) -> np.array:

    dims = []

    for frame in U.trajectory:

        dims.append(frame.dimensions[:3])

    return np.max(dims, axis = 0)

def CalcDensity(GROUP        : mda.AtomGroup,
                CENTER_GROUP : mda.AtomGroup,
                BINS         : np.array) -> np.array:


    out = np.zeros((BINS.shape[0] - 1))

    POS = GROUP.positions.copy()
    CEN = CENTER_GROUP.center_of_mass()

    POS[:,0] -= CEN[0]
    POS[:,1] -= CEN[1]
    POS[:,2] -= CEN[2]

    for i in range(len(BINS)-1):

        idx_a = POS[:,2] >= BINS[i]
        idx_b = POS[:,2] <  BINS[i+1]

        idx = idx_a * idx_b

        dens = np.sum(GROUP[idx].masses)/(DIMS[0]*DIMS[1]*(BINS[1]-BINS[0]))
        out[i] = dens

    return out


if __name__ == '__main__':

    INDEX_FILE = ReadNDX(options.ndx_file)

    U = mda.Universe(options.tpr_file,
                     options.xtc_file,
                     in_memory = False)
    
    
    if not os.path.exists(f'{options.path}/dims.pk'):
        DIMS = CalcBoxSize(U)
        with open(f'{options.path}/dims.pk', 'wb') as f:
            pk.dump(DIMS,f,
                    protocol = 4)
        f.close()
    else:
        with open(f'{options.path}/dims.pk', 'rb') as f:
            DIMS = pk.load(f)
        f.close()
        

    GROUP        = sum([U.select_atoms(f'index {idx}') for idx in INDEX_FILE[options.selection]])
    CENTER_GROUP = sum([U.select_atoms(f'index {idx}') for idx in INDEX_FILE[options.center]])

    BINS = np.arange(-DIMS[2]/2, DIMS[2]/2, 0.5)

    DENSITY = []

    t0 = U.trajectory[0].time

    idx_0 = int((options.initial_time - t0)/U.trajectory.dt)
    idx_1 = int((options.final_time - t0)/U.trajectory.dt)

    final_time = U.trajectory[-1].time

    for n, frame in enumerate(U.trajectory[idx_0:idx_1]):

        print(f'{int(100 * (n/(idx_1-idx_0-1)))}',
              end = '\r')

        d = CalcDensity(GROUP, CENTER_GROUP, BINS)
        DENSITY.append(d)

    DENSITY = np.array(DENSITY)

    X = 0.5*(BINS[1:] + BINS[:-1])

    DENSITY = Density(X, DENSITY)

    with open(f'{options.path}/{options.outfile}',
              'wb') as f:
        pk.dump(DENSITY, f,
                protocol = 4)
    f.close()















