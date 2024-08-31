#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 09:56:00 2022

@author: alexandre
"""

import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder as lf
import numpy as np
import pickle
from optparse import OptionParser

parser = OptionParser(usage       = 'python flipflop.py -p input.pdb -x input.xtc -r <resolution> -k Pickles -o output -n <n_proc>',
                      prog        = 'FlipFlop',
                      description = 'This program calculates the translocation'
                                    ' of lipids from one leaflet to the other.')

parser.add_option('-p', '--pdb',
                  action='store', type = 'string',
                  help = 'Path to input structure file')

parser.add_option('-x', '--xtc',
                  action='store', type = 'string',
                  help = 'Path to input trajectory file')

parser.add_option('-r', '--res',
                  action='store', type = 'string',
                  help = 'Resolution of the simulation (AA or CG)')

parser.add_option('-k', '--kle',
                  action='store', type = 'string',
                  help = 'Path to Pickles folder')

parser.add_option('-o', '--out',
                  action='store', type = 'string',
                  help = 'Output pickled file')

parser.add_option('-n', '--prc',
                  action='store', type = 'int',
                  help = 'Number of Processors, defaults to 6',
                  default = 6)

options, arguments = parser.parse_args()

def getDist(gx, gy, COM, P):
    '''
    Calculate the distance of the phosphate atoms to the membranes Center of 
    Mass as a measure of Membrane's thickess.

    Parameters
    ----------
    gx : NumPy array
        Grid spanning the X dimension.
    gy : NumPy array
        Grid spanning the Y dimension.
    COM : NumPy array
        2D array determining the center of the membrane at each (x,y) point.
    P : NumPy array
        Array containing the Phosphate positions.
    Z : NumPy array
        Output array to store distance of Phosphate to COM.
    l : int
        Integrer describing the analyzed leaflet (0 --> upper, 1--> lower).

    Returns
    -------
    None.

    '''

    dist = np.zeros((len(P)))

    for i, p in enumerate(P):
        
        setx = gx < p[0] # Gridx smaller than position X
        sety = gy < p[1] # Gridy smaller than position Y
        
        idx = np.where(setx==False)[0][0] - 1 # Get grid indices where it is 
        idy = np.where(sety==False)[0][0] - 1 # no longer smaller than (X,Y)
        
        com = COM[idx,idy] # Read center of mass at previous indices
        
        dist[i] += p[2]-com # Compute Z distance to COM
    
    return dist

#################### LOAD AUXILIARY TOPOGRAPHY AND GRID ARRAYS ##############

with open('topog.pickle', 'rb') as f:
    H = pickle.load(f)
f.close()

with open('gridX.pickle', 'rb') as f:
    gX = pickle.load(f)
f.close()

with open('gridY.pickle', 'rb') as f:
    gY = pickle.load(f)
f.close()

#############################################################################


################### LOAD GROMACS TRAJECTORY AND FIND LEAFLETS ###############

u = mda.Universe(options.pdb, options.xtc, in_memory = False)

res = {'AA':'P', 'CG':'PO4'}

P = u.select_atoms('name %s' % res[options.res], periodic = True, updating = True)

Leafs = lf(u, P, pbc = True)

Lf_up = Leafs.groups(0)
Lf_dn = Leafs.groups(1)

#############################################################################


################## SEARCH LOOP ##############################################

Z = np.zeros((u.trajectory.n_frames, 2, len(Lf_up)))

def Compute ( frame, H, gX, gY, Leaflet_up, Leaflet_down):

    h      = H[frame]
    gx, gy = gX[frame], gY[frame]

    COM = np.mean(h, axis = 0)

    Leaflet_up.universe.trajectory[frame]
    Leaflet_down.universe.trajectory[frame]

    Posup = Leaflet_up.positions
    Posdn = Leaflet_down.positions

    dist_up = getDist(gx, gy, COM, Posup)
    dist_dn = getDist(gx, gy, COM, Posdn)

    return dist_up, dist_dn


#############################################################################

if __name__=='__main__':

    import multiprocessing
    from functools import partial

    run_frames = partial(Compute, H = H, gX = gX, gY = gY, Leaflet_up = Lf_up, Leaflet_down = Lf_dn)

    Frames = np.arange(u.trajectory.n_frames)

    pool = multiprocessing.Pool(options.prc)

    Result = pool.map(run_frames, Frames)

    Z[:,0] += np.array([frame[0] for frame in Result])
    Z[:,1] += np.array([frame[1] for frame in Result])

################# OUTPUT ####################################################

with open('%s.pickle' % options.out, 'wb') as f:
    pickle.dump(Z,f, protocol = 4)
f.close()
