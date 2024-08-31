# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 12:28:24 2023

@author: alexandre
"""

import os
import numpy      as np
import pickle     as pk 
import MDAnalysis as mda

from scipy.integrate import simpson
from natsort         import natsorted as ns

global WATER_MASS
WATER_MASS = 2.99E-26 # kg / molecule

global strformat
strformat = '{:>3.2e}'

class Density:

    def __init__(self, X, DENSITY):

        self.X       = X
        self.DENSITY = DENSITY

def GetBoxArea(U : mda.Universe) -> np.array:

    OUT = []

    for frame in U.trajectory:

        area = (frame.dimensions[0] * frame.dimensions[1]) # Area in Angstrom^2

        OUT.append([frame.time, area])

    return np.array(OUT)

def Density2Particles(Density : np.array,
                      Area    : np.array,
                      MinTime : float,
                      MaxTime : float) -> np.array:

    idx_min = Area[Area[:,0] < MinTime].shape[0]
    idx_max = Area[Area[:,0] < MaxTime].shape[0]

    mean_area = Area[idx_min:idx_max,1].mean()

    Density[:,1] *= 0.05547068932594996 * mean_area # Final units are particles per Angstrom

def GetDensMaxima(Density : np.array) -> tuple:

    segment_1 = Density[:int(Density.shape[0]/2)]
    segment_2 = Density[int(Density.shape[0]/2):]

    idx_1 = np.where(segment_1[:,1] == segment_1[:,1].max())[0][0]
    idx_2 = np.where(segment_2[:,1] == segment_2[:,1].max())[0][0]

    return segment_1[idx_1,0], segment_2[idx_2,0]

def IntegrateInRange(Density : np.array,
                     Min_val : float,
                     Max_val : float) -> float:

    idx_min = Density[Density[:,0] <= Min_val].shape[0]
    idx_max = Density[Density[:,0] <= Max_val].shape[0]

    integral = simpson(Density[idx_min:idx_max,1], Density[idx_min:idx_max,0])

    return integral


if __name__ == '__main__':

    FILES_SOL = ns([_ for _ in os.listdir('./DensBootstrap/') if 'sol' in _])
    FILES_GLC = ns([_ for _ in os.listdir('./DensBootstrap/') if 'glc' in _])

    U = mda.Universe('prod.tpr',
                     '400_500ns.xtc')

    Area = GetBoxArea(U)

    DATA = []

    for dens_sol, dens_glc in zip(FILES_SOL, FILES_GLC):

        with open(f'./DensBootstrap/{dens_sol}', 'rb') as f:
            data_sol = pk.load(f)
            data_sol = np.array([data_sol.X, data_sol.DENSITY.mean(axis=0)]).T
        f.close()

        with open(f'./DensBootstrap/{dens_glc}', 'rb') as f:
            data_glc = pk.load(f)
            data_glc = np.array([data_glc.X, data_glc.DENSITY.mean(axis=0)]).T
        f.close()

        times = [float(_) for _ in dens_sol[-16:-3].split('_')]

        Density2Particles(data_sol, Area, times[0], times[1])

        MinZ, MaxZ = GetDensMaxima(data_glc)

        N_Waters = IntegrateInRange(data_sol, MinZ, MaxZ)

        DATA.append(N_Waters)

    print(f'Number of waters: {strformat.format(np.mean(DATA))} +- {strformat.format(np.std(DATA))}')
