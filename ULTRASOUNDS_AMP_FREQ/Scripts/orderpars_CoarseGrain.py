#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 17:21:50 2022

@author: alexandre
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder as lf
import pickle
from optparse import OptionParser

parser = OptionParser(usage       = 'python orderpars_CoarseGrain.py -p input.pdb -x input.xtc -i input.itp -k Pickles -o output -n <n_proc>',
                      prog        = 'OrderParameters CoarseGrain',
                      description = 'This programs calculates the CC (Segmental) Order'
                                    ' Parameters from simulations of Coarse Grained (Martini)'
                                    ' membrane bilayers. A structre (gro/pdb),'
                                    ' trajectory (trr/xtc), and topology (itp)'
                                    ' file must be provided. The program'
                                    ' assumes the user has performed the surface'
                                    ' determination step, and knows the path to the'
                                    ' corresponding Pickle directory,'
                                    ' otherwise it will not work.')

parser.add_option('-p', '--pdb',
                  action='store', type = 'string',
                  help = 'Path to input structure file')

parser.add_option('-x', '--xtc',
                  action='store', type = 'string',
                  help = 'Path to input trajectory file')

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

def unit_vector(vector):
    '''
    Makes a vector unit length

    Parameters
    ----------
    vector : NumPy array
        Vector to make unitary.

    Returns
    -------
    UnitaryVector
        Vector of unit length.

    '''

    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    '''
    Compute the angle between two n-dimensional vectors.

    Parameters
    ----------
    v1 : NumPy array
        First vector.
    v2 : NumPy array
        Second vector.

    Returns
    -------
    theta: float
        Angle between vectors.

    '''
    
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def OrderParameter(Phosphates, TailBeads, gx, gy, lrs):
    '''
    Compute the (not yet) Segmental Order Parameters of a lipid tail of interest.

    Parameters
    ----------
    Phosphates : MDAnalysis.AtomGroup
        Atom group containing the phosphates of the lipid of interest in one
        leaflet
    TailBeads : MDAnalysis.AtomGroup
        Atom group containing the beads of the desired tail in one 
        leaflet
    lrs : NumPy array
        Multidimensional array containing the Local Reference System at each
        point in the membrane surface

    Returns
    -------
    ORD : NumPy array
        2D array containing the (not yet) segmental order parameter of each 
        bead of the lipid tail

    '''
    startT = 0
    endT   = 0
    ORD = []
    for P in Phosphates:
        pos = P.position
        res = P.resid
        setx = gx < pos[0] # Gridx smaller than position X
        sety = gy < pos[1] # Gridy smaller than position Y
        idx = np.where(setx==False)[0][0]-1 # Get indices of grid where it
        idy = np.where(sety==False)[0][0]-1 # is no longer smaller than pX & pY
        n = lrs[idx,idy,2] # Normal vector at selected grid point
        for T in TailBeads[startT:]:
            resT = T.resid
            if resT == res:
                endT += 1
            else:
                break
        tail_atoms = TailBeads[startT:endT].positions
        vecs = tail_atoms[1:] - tail_atoms[:-1]
        Ordp = np.zeros((vecs.shape[0]))
        for v, vec in enumerate(vecs):
            theta   = angle_between(vec, n)
            ordp    = np.cos(theta)**2
            Ordp[v] = ordp
        ORD.append(Ordp)
        startT = endT
    ORD = np.array(ORD)
    return ORD

#################### LOAD AUXILIARY LRS AND GRID ARRAYS #####################

with open('./%s/lrss.pickle' % options.kle, 'rb') as f:
    LRS = pickle.load(f)
f.close()

with open('./%s/gridX.pickle' % options.kle, 'rb') as f:
    grX = pickle.load(f)
f.close()

with open('./%s/gridY.pickle' % options.kle, 'rb') as f:
    grY = pickle.load(f)
f.close()

#############################################################################


############ LOAD TRAJECTORY AND OBTAIN ATOMS OF INTEREST ##################

u = mda.Universe(options.pdb, options.xtc, in_memory = True)

notlipid = ['W', 'ION']

resnames = []

for res in u.atoms.residues:
    if not res.resname in resnames and not res.resname in notlipid:
        resnames.append(res.resname)

PO4   = u.select_atoms('name PO4', updating = True, periodic = True)

Leafs = lf(u, PO4, pbc = True)

Leaf_up = Leafs.groups(0)
Leaf_dn = Leafs.groups(1)

Phosphates = []

for residue in resnames:

    P_up = mda.AtomGroup([atom for atom in Leaf_up.atoms if atom.resname == residue])
    P_dn = mda.AtomGroup([atom for atom in Leaf_dn.atoms if atom.resname == residue])

    Phosphates.append([P_up,P_dn])

Resids_up = Leaf_up.residues
Resids_dn = Leaf_dn.residues

Resup = [res.resid for res in Resids_up]
Resdn = [res.resid for res in Resids_dn]

Tail1 = u.select_atoms('name C1A or name D2A or name C3A or name C4A', updating = True, periodic = True)
Tail2 = u.select_atoms('name C1B or name C2B or name C3B or name C4B', updating = True, periodic = True)

Tails = []

for residue in resnames:

    T1_up = mda.AtomGroup([atom for atom in Tail1.atoms if atom.resid in Resup and atom.resname == residue])
    T1_dn = mda.AtomGroup([atom for atom in Tail1.atoms if atom.resid in Resdn and atom.resname == residue])

    T2_up = mda.AtomGroup([atom for atom in Tail2.atoms if atom.resid in Resup and atom.resname == residue])
    T2_dn = mda.AtomGroup([atom for atom in Tail2.atoms if atom.resid in Resdn and atom.resname == residue])

    Tails.append([[T1_up,T1_dn],[T2_up,T2_dn]])



################### SEARCH LOOP #############################################
'''
ORDER = []

lt = len(u.trajectory)

for f, frame in enumerate(u.trajectory):
    
    print('Analyzed %i %% Frames' %(100*(f+1)/lt), end = '\r')
    
    gx, gy = grX[f], grY[f]
    lrs    = LRS[f]
    
    order = []

    for res in range(len(resnames)):

        OrdPars1_up = OrderParameter(Phosphates[res][0], Tails[res][0][0], gx, gy, lrs[0])
        OrdPars1_dn = OrderParameter(Phosphates[res][1], Tails[res][0][1], gx, gy, lrs[1])

        OrdPars2_up = OrderParameter(Phosphates[res][0], Tails[res][1][0], gx, gy, lrs[0])
        OrdPars2_dn = OrderParameter(Phosphates[res][1], Tails[res][1][1], gx, gy, lrs[1])

        order.append([[OrdPars1_up, OrdPars1_dn],[OrdPars2_up, OrdPars2_dn]])

    ORDER.append(order) 
    
ORDER = np.array(ORDER)
'''

ORDER = np.zeros((LRS.shape[0], len(resnames), 2, 2, len(Phosphates[0][0]), int(len(Tails[0][0][0])/len(Phosphates[0][0])-1)))

def Compute(frame, grX, grY, LRS):
    
    gx, gy = grX[frame], grY[frame]
    lrs    = LRS[frame]
    
    order = []

    for res in range(len(resnames)):

        Phosphates[res][0].universe.trajectory[frame]
        Phosphates[res][1].universe.trajectory[frame]

        Tails[res][0][0].universe.trajectory[frame]
        Tails[res][0][1].universe.trajectory[frame]
        Tails[res][1][0].universe.trajectory[frame]
        Tails[res][1][1].universe.trajectory[frame]

        OrdPars1_up = OrderParameter(Phosphates[res][0], Tails[res][0][0], gx, gy, lrs[0])
        OrdPars1_dn = OrderParameter(Phosphates[res][1], Tails[res][0][1], gx, gy, lrs[1])

        OrdPars2_up = OrderParameter(Phosphates[res][0], Tails[res][1][0], gx, gy, lrs[0])
        OrdPars2_dn = OrderParameter(Phosphates[res][1], Tails[res][1][1], gx, gy, lrs[1])

        order.append([[OrdPars1_up, OrdPars1_dn],[OrdPars2_up, OrdPars2_dn]])

    order = np.array(order)

    return order
 
#############################################################################

if __name__=='__main__':

    import multiprocessing
    from functools import partial

    run_frames = partial(Compute, grX = grX, grY = grY, LRS = LRS)

    Frames = np.arange(u.trajectory.n_frames)

    pool = multiprocessing.Pool(options.prc)

    Result = pool.map(run_frames, Frames)

    ORDER += np.array(Result)

############### OUTPUT ######################################################

for r, residue in enumerate(resnames):

    with open('%s_%s.pickle' % (options.out, residue), 'wb') as f:
        pickle.dump(ORDER[:,r], f, protocol = 4)
    f.close()
