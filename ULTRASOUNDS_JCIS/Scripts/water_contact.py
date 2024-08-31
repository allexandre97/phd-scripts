#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 12:31:43 2022

@author: alexandre
"""
import os
import MDAnalysis as mda
import numpy      as np
from   math       import sqrt
from   numba      import guvectorize
import pickle

@guvectorize(['float32[:,:], float32[:,:], float32, float32, float32, float32[:,:]'],
             '(n,m), (o,m), (), (), () -> (n,o)', target = 'parallel')
def DistanceMatrix(W_Pos, C_Pos, lx, ly, lz, DMat):
    
    for nw in range(W_Pos.shape[0]):
        for nc in range(C_Pos.shape[0]):
            
            dx = (W_Pos[nw,0] - C_Pos[nc,0]) # % lx
            dy = (W_Pos[nw,1] - C_Pos[nc,1]) # % ly
            dz = (W_Pos[nw,2] - C_Pos[nc,2]) # % lz
            DMat[nw,nc] = sqrt(dx**2 + dy**2 + dz**2)

def Compute(frame, Water, Tails):
    
    dims = Water.universe.trajectory[frame].dimensions
    lx, ly, lz = dims[0], dims[1], dims[2]
    
    Water.universe.trajectory[frame]
    
    W_Pos = Water.positions
     
    N = []
    
    DMat = np.zeros((W_Pos.shape[0], Tails[0][0].positions.shape[0]))
    
    for bead in range(4):
        
        Tails[0][bead].universe.trajectory[frame]
        Tails[1][bead].universe.trajectory[frame]
        
        C_A_Pos = Tails[0][bead].positions
        C_B_Pos = Tails[1][bead].positions
        
        DistanceMatrix(W_Pos, C_A_Pos, lx, ly, lz, DMat)
        n_contacts_A = np.where(DMat < 5)[0].shape[0]
                 
        DistanceMatrix(W_Pos, C_B_Pos, lx, ly, lz, DMat)
        n_contacts_B = np.where(DMat < 5)[0].shape[0]
                 
        N.append([n_contacts_A, n_contacts_B])

    return np.array(N)

        
u = mda.Universe('GRO/initial.gro', 'XTC/center.xtc', in_memory = True)

W = u.select_atoms('resname W', periodic = True, updating = True)

Tail_A = [u.select_atoms(f'name {t}{n}A', periodic = True, updating = True) for n, t in zip(range(1, 5), ['C', 'D', 'C', 'C'])]
Tail_B = [u.select_atoms(f'name C{n}B', periodic = True, updating = True) for n in range(1, 5)]

DMat = np.zeros((W.atoms.positions.shape[0], Tail_A[0].atoms.positions.shape[0]), dtype = 'float32')

NC = []

if __name__ == '__main__':
    
    import multiprocessing as mp
    from   functools       import partial
    
    chunksize = int(os.cpu_count()*10)
    MaxIter   = int(np.ceil(u.trajectory.n_frames/chunksize))
    
    with mp.Pool(os.cpu_count()) as pool:
        
        for Iter in range(MaxIter):
            beg = Iter*chunksize
            if (Iter+1)*chunksize < u.trajectory.n_frames:
                end = (Iter+1)*chunksize
            else:
                end = u.trajectory.n_frames
                
            run_frames = partial(Compute,
                                 Water    = W,
                                 Tails    = [Tail_A, Tail_B])
            
            Frames = np.arange(beg, end)

            result = pool.map(run_frames, Frames)
            
            for res in result:
                
                NC.append(res)
            
            print(f'Done {int(100*((Iter+1)/MaxIter))} %', end = '\r')
            
    pool.close()
    pool.join()

NC = np.array(NC)

with open('water_contact.pickle', 'wb') as f:
    
    pickle.dump(NC, f)

f.close()
