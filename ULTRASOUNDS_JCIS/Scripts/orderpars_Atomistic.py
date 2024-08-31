#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 14:47:05 2022

@author: alexandre
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder as lf
import fileinput
import pickle
from optparse import OptionParser

parser = OptionParser(usage       = 'python orderpars_Atomistic.py -p input.pdb -x input.xtc -i input.itp -k Pickles -o output -n <n_proc>',
                      prog        = 'OrderParameters Atomistic',
                      description = 'This programs calculates the CH Order'
                                    ' Parameters from simulations of atomistic'
                                    ' membrane bilayers. A structre (gro/pdb),'
                                    ' trajectory (trr/xtc), and topology (itp)'
                                    ' file must be provided. The program'
                                    ' assumes the user has performed the surface'
                                    ' determination step, and has the corresponding'
                                    ' pickle files in the Pickle directory,'
                                    ' otherwise it will not work')

parser.add_option('-p', '--pdb',
                  action='store', type = 'string',
                  help = 'Path to input structure file')

parser.add_option('-x', '--xtc',
                  action='store', type = 'string',
                  help = 'Path to input trajectory file')

parser.add_option('-i', '--itp',
                  action='store', type = 'string',
                  help = 'Path to input topology file')

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

def ReadITPfile(itp_file):
    '''
    Read a GROMACS topology file to obtain the atoms and their respective
    connectivities.

    Parameters
    ----------
    itp_file : String
        String containing the path to the desired topology file.

    Returns
    -------
    list
        List containing the atoms and the bonds between them.

    '''
    control=-1; atoms=[]; bonds=[]
    for line in fileinput.input(itp_file):
        if line.strip():
            fields=line.split()
            if (fields[0]=='[') and (fields[2]==']') and (fields[1]=="atomtypes" or fields[1]=="defaults"):
                control=0
            elif (fields[0]=='[') and (fields[2]==']') and (fields[1]=="atoms") :
                control=1
            elif (fields[0]=='[') and (fields[2]==']') and (fields[1]=="bonds") :
                control=2
            elif (fields[0]=='[') and (fields[2]==']') and (fields[1]=="pairs") :
                control=3
            elif (fields[0]=='[') and (fields[2]==']') and (fields[1]=="angles") :
                control=4
    	    # write to "atoms" and "bonds"
            if (control==1) and (fields[0][0]!=';') and (fields[0]!='[') and (len(fields)>1):
                atoms.append(fields)
            elif (control==2) and (fields[0][0]!=';') and (fields[0]!='[') and (len(fields)>1):
                bonds.append(fields)
    return [atoms, bonds]

def GetTailsSelections(itpfile, universe, leaflet, Tail):
    '''
    Return the ordered selection of carbons and hydrogens in a lipid tail of
    interest to perform order parameter analysis.
    
    Parameters
    ----------
    itpfile : string
        Topology files of the lipid to calculate order parameters for
    universe : MDAnalysis.Universe
        MDAnalysis universe object from which atoms will be selected.
    leaflet : MDAnalysis.AtomGroup.Residues
        Residue list all belonging to one specific leaflet
    Tail : string
        Name (SN1 or SN2) of the tail to study  

    Returns
    -------
    Carbons_leaflet : MDAnalysis.AtomGroup
        Atom group contining the carbon atoms of the tail to study
    Hydrogens_leaflet : MDAnalysus.AtomGroup
        Atom group containing the hydrogen atoms of the tail to study

    '''
    
    tailtype = {'SN1':2, 'SN2':3} # This dict translate Tail input to wht appears
                                  # in the topology
    
    Atoms_Bonds = ReadITPfile(itpfile)
    
    '''First we read the carbon atoms in the tail'''
    tail = []
    for atom in Atoms_Bonds[0]:
        if atom[4][0:2] == 'C%i' % tailtype[Tail] and len(atom[4]) > 2:
            tail.append([atom[0], atom[4]])
    
    '''Get their connectivities'''
    Connectivity = [] 
    for atom_tail in tail[1:]:
        idx     = int(atom_tail[0])
        Connect = [atom_tail[1]]
        for atom in Atoms_Bonds[1]:
            if int(atom[0]) == idx:
                connect = int(atom[1])
                Connect.append(Atoms_Bonds[0][connect-1][4])
        Connectivity.append(Connect)
    
    
    '''Generate selection strings for MDAnalysis'''
    carbons   = ''
    hydrogens = ''
    for c, carbon in enumerate(Connectivity):
        carbons   += 'name %s' % carbon[0]
        hydrogens += 'name %s' % carbon[1]
        if c < len(Connectivity)-1:
            carbons   += ' or '
            hydrogens += ' or '
    
    '''First select all atoms within selection string'''
    Carbons_all = universe.select_atoms(carbons, updating = True, periodic = True)
    Hydrogens_all = universe.select_atoms(hydrogens, updating = True, periodic = True)

    '''Filter first selection based on if their residue is present in the desired leaflet'''
    Carbons_leaflet   = mda.AtomGroup([atom for atom in Carbons_all if atom.resid in leaflet])
    Hydrogens_leaflet = mda.AtomGroup([atom for atom in Hydrogens_all if atom.resid in leaflet])
    
    return Carbons_leaflet, Hydrogens_leaflet

def GetOrderParameters(Phosphates, TailCarbons, TailHydrogens, gx, gy, lrs):
    '''
    Compute the (not yet) Order Parameters of a lipid tail of interest.

    Parameters
    ----------
    Phosphates : MDAnalysis.AtomGroup
        Atom group containing the phosphates of the lipid of interest in one
        leaflet
    TailCarbons : MDAnalysis.AtomGroup
        Atom group containing the carbon atoms of the desired tail in one 
        leaflet
    TailHydrogens : MDAnalysis.AtomGroup
        Atom group containing the hydrogen atoms of the desired tail in one
        leaflet
    lrs : NumPy array
        Multidimensional array containing the Local Reference System at each
        point in the membrane surface

    Returns
    -------
    ORD : NumPy array
        2D array containing the (not yet) order parameter of each carbon of
        the lipid tail

    '''
        
    startT, endT = 0, 0 # Auxiliary ints to delimitate tails
    ORD = []            # List to store order parameters
    for P in Phosphates:
        pos = P.position
        res = P.resid
        setx = gx < pos[0] # Gridx smaller than position X
        sety = gy < pos[1] # Gridy smaller than position Y
        
        idx = np.where(setx==False)[0][0]-1 # Get indices of grid where it
        idy = np.where(sety==False)[0][0]-1 # is no longer smaller than pX & pY
        
        n = lrs[idx,idy,2] # Normal vector at selected grid point
        
        for T in TailCarbons[startT:]: # Filter carbon atoms based on their
                                       # belonging to same residue than the
                                       # phosphate used to determine position
            resT = T.resid
            if resT == res:
                endT += 1
            else:
                break
            
        tail_carbons   = TailCarbons[startT:endT].positions   # Select pertinent atoms
        tail_hydrogens = TailHydrogens[startT:endT].positions
        vecs = tail_hydrogens - tail_carbons  # Calculate vectors
        
        ordp = np.zeros((vecs.shape[0])) # Initialize array to store ordp.
        
        '''Calculate ordp for each vector'''
        for v, vec in enumerate(vecs):
            theta    = angle_between(vec, n)
            ordp[v]  = np.cos(theta)**2
            
        ORD.append(ordp)
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

############################################################################


############ LOAD TRAJECTORY AND OBTAIN ATOMS OF INTEREST ##################

residue = options.itp[-8:-4]

u = mda.Universe(options.pdb, options.xtc, in_memory = True)

P = u.select_atoms('name P and resname %s' % residue, periodic = True, updating = True)

Leaflets = lf(u, P, pbc = True)

Leaflet_up = Leaflets.groups(0)
Leaflet_dn = Leaflets.groups(1)

Residues_up = Leaflet_up.residues
Residues_dn = Leaflet_dn.residues

Resup = [res.resid for res in Residues_up]
Resdn = [res.resid for res in Residues_dn]

Carbons1_up, Hydrogens1_up = GetTailsSelections(options.itp, u, Resup, 'SN1')
Carbons1_dn, Hydrogens1_dn = GetTailsSelections(options.itp, u, Resdn, 'SN1')

Carbons2_up, Hydrogens2_up = GetTailsSelections(options.itp, u, Resup, 'SN2')
Carbons2_dn, Hydrogens2_dn = GetTailsSelections(options.itp, u, Resdn, 'SN2')

#############################################################################


################### SEARCH LOOP #############################################

ratio1 = int(len(Carbons1_up)/len(Leaflet_up))

ratio2 = int(len(Carbons2_up)/len(Leaflet_up))

ORDER_SN1 = np.zeros((LRS.shape[0], 2, len(Leaflet_up), ratio1))
ORDER_SN2 = np.zeros((LRS.shape[0], 2, len(Leaflet_up), ratio2))

def Compute(frame, LRS, grX, grY):
    
    Leaflet_up.universe.trajectory[frame]
    Leaflet_dn.universe.trajectory[frame]

    Carbons1_up.universe.trajectory[frame]
    Carbons2_up.universe.trajectory[frame]

    Carbons1_dn.universe.trajectory[frame]
    Carbons2_dn.universe.trajectory[frame]

    lrs = LRS[frame]
    gx  = grX[frame]
    gy  = grY[frame] 
    
    OrdPars1_up = GetOrderParameters(Leaflet_up, Carbons1_up, Hydrogens1_up, gx, gy, lrs[0])
    OrdPars1_dn = GetOrderParameters(Leaflet_dn, Carbons1_dn, Hydrogens1_dn, gx, gy, lrs[1])

    OrdPars2_up = GetOrderParameters(Leaflet_up, Carbons2_up, Hydrogens2_up, gx, gy, lrs[0])
    OrdPars2_dn = GetOrderParameters(Leaflet_dn, Carbons2_dn, Hydrogens2_dn, gx, gy, lrs[1])

    return OrdPars1_up, OrdPars1_dn, OrdPars2_up, OrdPars2_dn
  
#############################################################################

if __name__=='__main__':

    import multiprocessing
    from functools import partial

    run_frames = partial(Compute, LRS = LRS, grX = grX, grY = grY)

    Frames = np.arange(LRS.shape[0])

    pool = multiprocessing.Pool(options.prc)

    Result = pool.map(run_frames, Frames)

    ORDER_SN1[:,0] += np.array([frame[0] for frame in Result])
    ORDER_SN1[:,1] += np.array([frame[1] for frame in Result])

    ORDER_SN2[:,0] += np.array([frame[2] for frame in Result])
    ORDER_SN2[:,1] += np.array([frame[3] for frame in Result])

############### OUTPUT ######################################################

with open(options.out+'_s1_%s.pickle' % residue, 'wb') as f:
    pickle.dump(ORDER_SN1, f)
f.close()

with open(options.out+'_sn2_%s.pickle' % residue, 'wb') as f:
    pickle.dump(ORDER_SN2, f)
f.close()

