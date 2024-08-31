# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 22:05:09 2023

@author: alexandre
"""

import os
import numpy             as np
import pickle            as pk
import MDAnalysis        as mda
#import matplotlib.pyplot as plt


from numba import guvectorize
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-n', '--ndx_file',
                  default = 'NDX/index_leaflet.ndx',
                  type    = str,
                  help    = 'input ndx file',
                  action  = 'store')
parser.add_option('-t', '--tpr_file',
                  default = 'prod.tpr',
                  type    = str,
                  help    = 'input tpr file',
                  action  = 'store')
parser.add_option('-x', '--xtc_file',
                  default = 'XTC/400_500.xtc',
                  type    = str,
                  help    = 'input xtc file',
                  action  = 'store')
parser.add_option('-s', '--selection',
                  default = 'POPE_UP',
                  type    = str,
                  help    = 'selection group',
                  action  = 'store')
parser.add_option('-i', '--initial_time',
                  default = 400000,
                  type    = int,
                  help    = 'initial time',
                  action  = 'store')
parser.add_option('-f', '--final_time',
                  default = 500000,
                  type    = int,
                  help    = 'final time',
                  action  = 'store')
parser.add_option('-p', '--path',
                  default = './POPE-POPG_3-1/Charmm/',
                  type    = str,
                  help    = 'master path',
                  action  = 'store')
parser.add_option('-o', '--outfile',
                  default = 'DNS/frames/frame',
                  type    = str,
                  help    = 'file name',
                  action  = 'store')
options, arguments = parser.parse_args()

class Density:
    
    def __init__(self, X, Y, DENSITY):
        
        self.X = X
        self.Y = Y
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

@guvectorize(['(float64[:], float64[:], float64[:,:], float64[:], float64[:,:])'],
              '(n),         (n),        (m,l),        (m)      -> (n,n)',
              target   = 'parallel',
              nopython = True,
              fastmath = True)
def CalcDens(X    : np.array, 
             Y    : np.array, 
             POS  : np.array,
             MASS : np.array,
             OUT  : np.array) -> np.array:
    
    A   = X[1] * Y[1]
        
    for i in range(X.shape[0] - 1):
        
        x_a = X[i]
        x_b = X[i+1]
        
        for j in range(Y.shape[0] - 1):
            
            y_a = Y[j]
            y_b = Y[j+1]
            
            for m, p in zip(MASS, POS):
                px = p[0]
                py = p[1]
                
                if x_a <= px and px < x_b and y_a <= py and py < y_b:
                    
                    OUT[i,j] += m 
            
            OUT[i,j] /= A
            


if __name__ == '__main__':
    
    INDEX_FILE = ReadNDX(f'{options.path}/{options.ndx_file}')
    

    U = mda.Universe(f'{options.path}/{options.tpr_file}',
                     f'{options.path}/{options.xtc_file}')
    
    GROUP = sum([U.select_atoms(f'index {idx}')
                 for idx in INDEX_FILE[options.selection]])
    
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
        
    X, Y = np.linspace(0.0, DIMS[0], 100), np.linspace(0.0, DIMS[1], 100)
    _ = 0.5 * (X[1:] + X[:-1])
    XX, YY = np.meshgrid(_, _)
    XX /= 10.
    YY /= 10.
    
    DENSITY = []

    t0 = U.trajectory[0].time

    idx_0 = int((options.initial_time - t0)/U.trajectory.dt)
    idx_1 = int((options.final_time - t0)/U.trajectory.dt)

    final_time = U.trajectory[-1].time
    d = np.zeros((X.shape[0],
                  Y.shape[0]))
    for n, frame in enumerate(U.trajectory[idx_0:idx_1]):

        print(f'{int(100 * (n/(idx_1-idx_0-1)))}',
              end = '\r')
        
        CalcDens(X, Y, GROUP.positions, GROUP.masses, d)
        
        DENSITY.append(d[:-1,:-1].copy())
        
        d[:,:] = 0.0
        
        
    print(f'DONE LOOP {options.path}')
    DENSITY = np.array(DENSITY)
    
    DENSITY = Density(X, Y, DENSITY.mean(axis = 0))
    
    with open(f'{options.path}/DNS/density_{options.selection}_{options.initial_time}.pk', 'wb') as f:
        pk.dump(DENSITY, f,
                protocol = 4)
    f.close()
    
    print(f'DONE WRITING {options.path}')