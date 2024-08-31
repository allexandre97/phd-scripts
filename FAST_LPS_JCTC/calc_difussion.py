# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 10:14:36 2023

@author: alexandre
"""


import numpy             as np
import pickle            as pk
import cblind            as cb
import MDAnalysis        as mda
import statsmodels.api   as sm
import matplotlib.pyplot as plt

from optparse            import OptionParser

class Velocities:
    
    def __init__(self,
                 BOUND,
                 BOUNDF,
                 VELS):
        
        self.BOUNDF = BOUNDF
        self.BOUND  = BOUND
        self.VELS   = VELS

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

parser.add_option('-o', '--outfile',
                  action = 'store', type = 'string',
                  default = 'diffusion.png',
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

LPS_N_ATOMS = 42
KB          = 8.314462618*(1E10/1E12)#1.380649E-23
T           = 313


def gen_COM_positions(AtomGroup, Atoms_Per_Molecule, N_Molecules):
    PCOM = []
    for n in range(N_Molecules):
        tmp = AtomGroup[n*Atoms_Per_Molecule:(n+1)*Atoms_Per_Molecule]
        PCOM.append(tmp.center_of_mass())
    PCOM = np.array(PCOM)
    return PCOM


def GetVels(AtomGroup,
            DIMS,
            frame,
            w,
            dt,
            Atoms_Per_Molecule,
            N_Molecules):
    
    AtomGroup.universe.trajectory[frame]
    POS0 = gen_COM_positions(AtomGroup, Atoms_Per_Molecule, N_Molecules)
    
    AtomGroup.universe.trajectory[frame+w]
    POS1 = gen_COM_positions(AtomGroup, Atoms_Per_Molecule, N_Molecules)
    
    dX = POS1 - POS0
    
    dX[:,0] -= np.round(dX[:,0]/DIMS[0])*DIMS[0]
    dX[:,1] -= np.round(dX[:,1]/DIMS[1])*DIMS[1]
    dX[:,2] -= np.round(dX[:,2]/DIMS[2])*DIMS[2]
    
    dXdt = dX/dt
    
    return dXdt
 

def GetBoundUnbound(BOUND, VELS):
    
    wherebound   = np.where(BOUND[:-1])
    whereunbound = np.where(np.logical_not(BOUND[:-1]))
    
    vels_bound   = VELS[wherebound[0], wherebound[1]]
    vels_unbound = VELS[whereunbound[0], whereunbound[1]]
    
    v3d_bound   = np.sqrt(vels_bound[:,0]**2 + vels_bound[:,1]**2 + vels_bound[:,2]**2)
    v3d_unbound = np.sqrt(vels_unbound[:,0]**2 + vels_unbound[:,1]**2 + vels_unbound[:,2]**2)
    
    return v3d_bound, v3d_unbound


if __name__ == '__main__':
    
    
    
    U = mda.Universe(options.tpr, options.xtc,
                     in_memory = False)
    
    FLPS_All     = U.select_atoms('resname FEMP',
                                  updating = True, periodic = True)
    
    FLPS_Bound   = U.select_atoms('byres (resname FEMP and around 6 resname PMB1)',
                                  periodic = True, updating = True)
    
    LPS_All     = U.select_atoms('resname REMP',
                                 updating = True, periodic = True)
    
    LPS_Bound   = U.select_atoms('byres (resname REMP and around 6 resname PMB1)',
                                 periodic = True, updating = True)
    
    
    MASS = np.sum(LPS_All[:LPS_N_ATOMS].masses)

    N_LPS  = int(len(LPS_All)/LPS_N_ATOMS)
    N_FLPS = int(len(FLPS_All)/LPS_N_ATOMS)
    BOUND = []
    BOUNDF = []
    
    dt = 20
    w  = 1#int(np.ceil(dt/U.trajectory.dt))

    print(w)

    idx_i =int(np.ceil(options.initial/U.trajectory.dt))
    idx_f =int(np.ceil(options.final/U.trajectory.dt))
    
    for frame in U.trajectory[idx_i:idx_f-w]:
        
        if N_LPS != 0:
        
            tmp = np.zeros((N_LPS), dtype = bool)
            ids = (np.unique(LPS_Bound.resids) - LPS_Bound.resids[0]) - 1
            tmp[ids] = True
            BOUND.append(tmp)
            
        if N_FLPS != 0:
            tmp = np.zeros((N_FLPS), dtype = bool)
            ids = (np.unique(FLPS_Bound.resids) - FLPS_Bound.resids[0]) - 1
            tmp[ids] = True
            BOUNDF.append(tmp)
    
    BOUND  = np.array(BOUND)
    BOUNDF = np.array(BOUNDF)
    
    VELS  = []

    for frame in range(idx_i, idx_f-w):
        
        print(f'Done: {int(100*(frame-idx_i)/((idx_f-idx_i)-w-1))} %', end = '\r')
        
        DIMS = U.trajectory[frame].dimensions[:3]
        
        if N_LPS != 0:
            VelsREMP = GetVels(LPS_All, DIMS, frame, w, dt, LPS_N_ATOMS, N_LPS)
        else:
            VelsREMP = np.zeros((N_FLPS, 3))
            
        if N_FLPS != 0:
            VelsFEMP = GetVels(FLPS_All, DIMS, frame, w, dt, LPS_N_ATOMS, N_FLPS)
        else:
            VelsFEMP = np.zeros((N_LPS, 3))
            
        VELS.append([VelsREMP, VelsFEMP])
    
    VELS = Velocities(BOUND, BOUNDF, np.array(VELS))
    
    with open('velocities.pk', 'wb') as f:
        pk.dump(VELS, f, protocol = 4)
    f.close()
    
