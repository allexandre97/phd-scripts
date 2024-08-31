# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 09:40:27 2023

@author: alexandre
"""

import numpy      as np
import pickle     as pk
import MDAnalysis as mda

from math     import sqrt
from numba    import cuda, guvectorize
from optparse import OptionParser

parser = OptionParser(usage       = 'python calc_clusters.py',
                      prog        = 'CalcClusters',
                      description = 'This is a program used to calculate clustering of molecules.') 

parser.add_option('-t', '--tpr', 
                  action = 'store', type = 'string',
                  default = './AllFast/5us.tpr',
                  help = 'Input .tpr File')

parser.add_option('-x', '--xtc',
                  action = 'store', type = 'string',
                  default = './AllFast/5us.xtc',
                  help = 'Input .xtc File')

parser.add_option('-o', '--out',
                  action = 'store', type = 'string',
                  default = 'clusters',
                  help = 'Output file master name')

parser.add_option('-s', '--start',
                  action = 'store', type = float,
                  default = 0.0,
                  help = 'Start time in ns')

parser.add_option('-f', '--finish',
                  action = 'store', type = float,
                  default = 300,
                  help = 'Finish time in ns')

parser.add_option('-a', '--grpa',
                  action = 'store', type = 'string',
                  default = 'resname PMB1',
                  help = 'Group A selection')

parser.add_option('-b', '--grpb',
                  action = 'store', type = 'string',
                  default = 'resname PMB1',
                  help = 'Group B selection')

options, args = parser.parse_args()

class Clusters:
    
    def __init__(self,
                 GroupA : str, GroupB : str,
                 Time : np.array,
                 ContactMatrix : np.array):
        
        self.GroupA        = GroupA
        self.GroupB        = GroupB
        self.Time          = Time
        self.ContactMatrix = ContactMatrix

def CountAtomsPerResidue(GROUP):
    N_Atoms = 0
    RES_0   = GROUP[0].resid
    for atom in GROUP:
        if atom.resid == RES_0:
            N_Atoms += 1
        else:
            break
    return N_Atoms

def BuildByUnit(Atoms, NAtoms):

    byunit = []
    tmp = [Atoms[0]]

    for atom in Atoms[1:]:
        if len(tmp) % NAtoms == 0:
            byunit.append(mda.AtomGroup(tmp))
            tmp = [atom]
        else:
            tmp.append(atom)
    byunit.append(tmp)

    return byunit

@cuda.jit(debug=False, opt=True)
def Dist(X, Y, L, D):

    width = X.shape[0]
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


@guvectorize(['(float32[:,:], float32, int32[:,:], int32[:,:])'],
             ' (l,n),         (),      (m,m)       -> (m,m)',
             target='parallel')
def Count(DistMat, Cutoff, dummy, N):
    
    X = DistMat.shape[0]
    Y = DistMat.shape[1]
    
    I = int(X/N.shape[0])
    J = int(Y/N.shape[1])
    
    i = 0
    for x in range(1, X+1):
        if ((x % I) == 0):
            i += 1
        j = 0
        for y in range(1, Y+1):
            if ((y % J) == 0):
                j += 1
            
            if DistMat[x-1,y-1] < Cutoff and (i != j):
                N[i-1,j-1] += 1
            
            

            

def Compute(GROUP_A, GROUP_B,
            DIM, Cutoff, blocks):
    
    D = np.zeros((len(GROUP_A), len(GROUP_B)),
                 dtype=np.float32)
    
    b   = blocks
    blc = (b, b)
    gx  = int(np.ceil(D.shape[0]/b))
    gy  = int(np.ceil(D.shape[1]/b))
    grd = (gx, gy)
    
    X_d   = cuda.to_device(GROUP_A.positions.copy())
    Y_d   = cuda.to_device(GROUP_B.positions.copy())
    D_d   = cuda.to_device(D)
    Dim_d = cuda.to_device(DIM)
    Dist[grd, blc](X_d, Y_d, Dim_d, D_d)
    D = D_d.copy_to_host()
    
    N_A = CountAtomsPerResidue(GROUP_A)
    N_B = CountAtomsPerResidue(GROUP_B)
    
    N = np.zeros((int(D.shape[0]/N_A), int(D.shape[1]/N_B)),
                 dtype = np.int32)
    
    Count(D, Cutoff, N, N)
    
    return N

if __name__ == '__main__':
    
    U = mda.Universe(options.tpr, options.xtc,
                     in_memory = True)
    
    GROUP_A = U.select_atoms(options.grpa,
                             periodic = True,
                             updating = True)
    
    GROUP_B = U.select_atoms(options.grpb,
                             periodic = True,
                             updating = True)
    
    RESNAME_A = GROUP_A[0].resname
    RESNAME_B = GROUP_B[0].resname
    
    
    timestep = U.trajectory.dt
    
    idx_start  = int(np.ceil( 1000 * options.start  / timestep ))
    idx_finish = int(np.ceil( 1000 * options.finish / timestep ))
    
    n_frames = idx_finish - idx_start - 1
    
    DATA = []
    
    for f, frame in enumerate(U.trajectory[idx_start:idx_finish+1]):
        
        N = Compute(GROUP_A, GROUP_B, frame.dimensions[:3], 6, 16)
        DATA.append(N)
        
        print(f'Done {int(100 * f/(n_frames))} %', end = '\r')
        
    
    DATA = np.array(DATA)
    
    T = np.arange(U.trajectory[idx_start].time,
                  U.trajectory[idx_finish].time + timestep,
                  timestep)
    
    ClusterData = Clusters(RESNAME_A, RESNAME_B, T, DATA)
    
    with open(f'{options.out}_{RESNAME_A}_{RESNAME_B}.pk', 'wb') as f:
        pk.dump(ClusterData, f,
                protocol = 4)
    f.close()