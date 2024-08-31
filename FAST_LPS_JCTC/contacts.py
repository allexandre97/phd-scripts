# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 12:02:01 2023

@author: alexandre
"""

import numpy             as np
import pickle            as pk
import MDAnalysis        as mda

from math                import sqrt
from numba               import cuda

from optparse            import OptionParser

class Contacts:

    def __init__(self,
                 RESNAME1,
                 RESNAME2,
                 GRP1_Names,
                 GRP2_Names,
                 CONTACTS):

        self.RESNAME1   = RESNAME1
        self.RESNAME2   = RESNAME2
        self.GRP1_Names = GRP1_Names
        self.GRP2_Names = GRP2_Names
        self.CONTACTS   = CONTACTS


parser = OptionParser(usage       = 'python wavesim.py -i <input build file> -o <output name>',
                      prog        = 'WaveSim',
                      description = 'This is a program used to simulate propagating waves through a medium.')

parser.add_option('-t', '--tpr',
                  action = 'store', type = 'string',
                  default = '5us.tpr',
                  help = 'Input .tpr File')

parser.add_option('-x', '--xtc',
                  action = 'store', type = 'string',
                  default = '5us.xtc',
                  help = 'Input .xtc File')

parser.add_option('-a', '--grp1',
                  action = 'store', type = 'string',
                  default = 'resname PMB1',
                  help = 'Group A')

parser.add_option('-b', '--grp2',
                  action = 'store', type = 'string',
                  default = 'resname FEMP',
                  help = 'Output .png file')

parser.add_option('-o', '--outfile',
                  action = 'store', type = 'string',
                  default = 'halffast',
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

@cuda.jit(fastmath = True)
def Dist(X, Y, L, D):

    width  = X.shape[0]
    height = Y.shape[0]

    startX, startY = cuda.grid(2)
    gridX = cuda.gridDim.x * cuda.blockDim.x
    gridY = cuda.gridDim.y * cuda.blockDim.y

    for i in range(startX, width, gridX):

        x = X[i]

        for j in range(startY, height, gridY):

            y = Y[j]

            dx = x[0]-y[0]
            dy = x[1]-y[1]
            dz = x[2]-y[2]

            dx -= round(dx/L[0])*L[0]
            dy -= round(dy/L[1])*L[1]
            dz -= round(dz/L[2])*L[2]

            D[i, j] = sqrt(dx*dx + dy*dy + dz*dz)

def Compute(X, Y, D, DIM, blocks = 8):

    blc = (blocks, blocks)
    gx  = int(np.ceil(D.shape[0]/blocks))
    gy  = int(np.ceil(D.shape[1]/blocks))
    grd = (gx, gy)

    X_d = cuda.to_device(X)
    Y_d = cuda.to_device(Y)
    DIM_d = cuda.to_device(DIM)
    D_d = cuda.to_device(D)
    Dist[grd, blc](X_d, Y_d, DIM_d, D_d)
    D = D_d.copy_to_host()

    return D

if __name__ == '__main__':

    U = mda.Universe(options.tpr, options.xtc,
                     in_memory = False)

    GRP1  = U.select_atoms(options.grp1,
                          periodic = True, updating = True)
    GRP2 = U.select_atoms(options.grp2,
                          periodic = True, updating = True)


    RESNAME1 = GRP1.resnames[0]
    RESNAME2 = GRP2.resnames[0]

    GRP1_Names = np.unique(GRP1.names)
    GRP2_Names = np.unique(GRP2.names)

    CONTACTS = []

    idx_i =int(np.ceil(options.initial/U.trajectory.dt))
    idx_f =int(np.ceil(options.final/U.trajectory.dt))

    NFRAMES = idx_f - idx_i

    for f, frame in enumerate(U.trajectory[idx_i:idx_f]):

        print(f'Done {int(100*f/(NFRAMES-1))} %', end = '\r')

        NameMatrix = np.zeros((len(GRP1_Names),
                               len(GRP2_Names)))

        D = np.zeros((len(GRP1),
                      len(GRP2)))

        D = Compute(GRP1.positions, GRP2.positions, D, U.dimensions)

        where_contact = np.where(D < 6)

        for pair in zip(where_contact[0], where_contact[1]):

            idx_grp1 = np.where(GRP1_Names == GRP1.names[pair[0]])[0][0]
            idx_grp2 = np.where(GRP2_Names == GRP2.names[pair[1]])[0][0]

            NameMatrix[idx_grp1, idx_grp2] += 1

        CONTACTS.append(NameMatrix)

    CONTACTS = Contacts(RESNAME1, RESNAME2, GRP1_Names, GRP2_Names, np.array(CONTACTS))

    with open(f'contacts_{RESNAME1}_{RESNAME2}.pk', 'wb') as f:
        pk.dump(CONTACTS, f, protocol = 4)
    f.close()