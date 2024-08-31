#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 14:33:40 2022

@author: alexandre
"""

import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder as lf
import numpy as np
import numpy.matlib as npm
from scipy.spatial.transform import Rotation as rot
from optparse import OptionParser
import pickle

parser = OptionParser(usage       = 'python memsearch.py -p input.pdb -x input.xtc -r <resolution> -n <n_proc>',
                      prog        = 'MemSearch',
                      description = 'This program performs a grid search over.'
                                    ' the membrane to determine it surface, as well'
                                    ' as obtain the positions and lipid directors, both in grid'
                                    ' format and in individual lipid format.')

parser.add_option('-p', '--pdb',
                  action='store', type = 'string',
                  help = 'Path to input structure file')

parser.add_option('-x', '--xtc',
                  action='store', type = 'string',
                  help = 'Path to input trajectory file')

parser.add_option('-r', '--res',
                  action='store', type = 'string',
                  help = 'Resolution of the simulation (AA or CG)')

parser.add_option('-n', '--prc',
                  action='store', type = 'int',
                  help = 'Number of Processors, defaults to 6',
                  default = 6)


options, arguments = parser.parse_args()

# Q is a Nx4 numpy matrix and contains the quaternions to average in the rows.
# The quaternions are arranged as (w,x,y,z), with w being the scalar
# The result will be the average quaternion of the input. Note that the signs
# of the output quaternion can be reversed, since q and -q describe the same orientation
def averageQuaternions(Q):
    '''
    Compute the average of a series of quaternions.

    Parameters
    ----------
    Q : NumPy array
        Array (N,4) containing N quaternions.

    Returns
    -------
    AvQ : NumPy array
        Array (1,4) containing the average quaternion.

    '''
    
    M = Q.shape[0] # Number of quaternions to average
    A = npm.zeros(shape=(4,4))
    for i in range(0,M):
        q = Q[i,:]
        # multiply q with its transposed version q' and add A
        A = np.outer(q,q) + A
    # scale
    A = (1.0/M)*A
    # compute eigenvalues and -vectors
    eigenValues, eigenVectors = np.linalg.eig(A)
    # Sort by largest eigenvalue
    eigenVectors = eigenVectors[:,eigenValues.argsort()[::-1]]
    # return the real part of the largest eigenvector (has only real part)
    return np.real(eigenVectors[:,0].A1)

def Director(P, TA, TB):
    '''
    Obtain the Lipid Director Matrix based on the orientation of the tails.
    
    Given the two lipid tails as two vector V1 and V2, the director matrix has
    dimensions 3x3 where the first row is:
    
    ---> A = (V1 + V2) / |V1 + V2|
    
    the second row:
        
    ---> B = (V1 x V2) / |V1 x V2|
    
    and the third row:
    
    ---> C = (A x B) / |A x B|

    Parameters
    ----------
    P : MDAnalysis.AtomGroup
        Atom group containing the Phosphates of the lipids.
    TA : MDAnalysis.AtomGroup
        Atom group containing the last bead of the first tail of the lipids.
    TB : MDAnalysis.AtomGroup
        Atom group containing the last bead of the second tail of the lipids.

    Returns
    -------
    D : NumPy array
        Array (N,3,3) containing the Director Matrix for the N lipids.

    '''

    V1 = TA.positions - P.positions # Calculate the tail vectors as the
    V2 = TB.positions - P.positions # position of the last beads minus the
                                    # position of the phosphates
                                    
    D = np.zeros((P.positions.shape[0], 3, 3)) # Initialize directo matrices
    
    
    # Calculate the pertinent vectors
    
    A = V1+V2
    l = np.sqrt(A[:,0]**2+A[:,1]**2+A[:,2]**2)
    A[:,0] /= l
    A[:,1] /= l
    A[:,2] /= l
    
    B = np.cross(V1, V2)
    l = np.sqrt(B[:,0]**2+B[:,1]**2+B[:,2]**2)
    B[:,0] /= l
    B[:,1] /= l
    B[:,2] /= l
    
    C = np.cross(A, B)
    l = np.sqrt(C[:,0]**2+C[:,1]**2+C[:,2]**2)
    C[:,0] /= l
    C[:,1] /= l
    C[:,2] /= l
    
    D[:,0] += A
    D[:,1] += B
    D[:,2] += C
    
    return D

def Search(n, P, Dr, gX, gY):
    '''
    Perform a search over a grid defined over the surface of the membrane and
    return the average director and elevation at each point.

    Parameters
    ----------
    n : Int
        Number of grid points in X and Y dimension.
    P : NumPy array
        Array (N,3) containing the phosphate positions for the N lipids.
    Dr : NumPy array
        Array (N,3,3) containing the lipid director matrix for the N lipids.
    gX : NumPy array
        Array containing the bin edges of the grid in X.
    gY : NumPy array
        Array containing the bin edges of the grid in Y.

    Returns
    -------
    H : NumPy array
        Array (n,n) containing the surface of the membrane.
    D : NumPy array
        Array (n,n,3,3) containing the average directors over the surface of
        the grid.

    '''
    
    H = np.zeros((n, n))        # Initialize output matrices
    D = np.zeros((n, n, 3, 3))
    for i in range(n):
        
        seta = gX[i] <= P[:,0]  # All the X points bigger or equal then the gX point
        setb = gX[i+1] > P[:,0] # All the X points smaller than the gX+1 point
        setx = seta*setb        # Intersection of previous sets
        
        for j in range(n):
            seta = gY[j] <= P[:,1]
            setb = gY[j+1] > P[:,1] # Same protocol carried for Y dimension
            sety = seta*setb
            
            SET = setx*sety         # Intersection of X and Y sets
            
            if np.any(SET) == True: # If total set is not empty
            
                h = np.mean(P[SET,2], axis = 0) # Elevation is mean of phosphates Z coord.
                r = rot.from_matrix(Dr[SET])    # Load director matrix as rotation
                q = r.as_quat()                 # and translate it to quaternion 
                
                Q = averageQuaternions(q)       # Average the quaternions
                r = rot.from_quat(Q)            # load them as rotation
                R = r.as_matrix()               # and back to matrix form
                
                H[i,j] += h
                D[i,j] += R
    return H, D    

def Interpolate(H, n):
    '''
    Subsitute empty bins with a linear interpolation of the contiguous cells.

    Parameters
    ----------
    H : NumPy array
        Array (X,Y) containing the elevation function over the surface of the
        membrane.
    n : Int
        Lateral dimension of the H grid.

    Returns
    -------
    None.

    '''

    nan = np.array(np.where(H==0.))
    
    for nor in range(nan.shape[1]):
        
        i = nan[0,nor] # Get index of empty value
        j = nan[1,nor]
        
        i1  = (i + 1) % n # Prevent overflow
        i_1 = i - 1
        j1  = (j + 1) % n
        j_1 = j - 1
        
        H[i,j] = (H[i_1,j] + H[i1,j] + H[i,j_1] + H[i,j1])/4 # Substitute
        

############### LOAD TRAJECTORY, LEAFLETS AND RELEVANT ATOMS ################    

u = mda.Universe(options.pdb, options.xtc, in_memory = False)

resP = {'AA':'P', 'CG':'PO4'}

P = u.select_atoms('name %s' % resP[options.res], periodic = True, updating = True)

leafs = lf(u, P, pbc = True)
Lf_up = leafs.groups(0)
Lf_dn = leafs.groups(1)

Res_up = Lf_up.residues
Res_dn = Lf_dn.residues

RN_up = []
RN_dn = []
for ru, rd in zip(Res_up, Res_dn):
    RN_up.append(ru.resname)
    RN_dn.append(rd.resname)
RN = np.array([RN_up, RN_dn])

with open('resnames.pickle', 'wb') as f:
    pickle.dump(RN,f)
f.close()

P_up = Lf_up.positions
P_dn = Lf_dn.positions

resT1 = {'AA':'C218', 'CG':'C4A'}
resT2 = {'AA':'C316', 'CG':'C4B'}

TA_up = mda.AtomGroup([atom for atom in Res_up.atoms if atom.name == resT1[options.res]])
TA_dn = mda.AtomGroup([atom for atom in Res_dn.atoms if atom.name == resT1[options.res]])
TB_up = mda.AtomGroup([atom for atom in Res_up.atoms if atom.name == resT2[options.res]])
TB_dn = mda.AtomGroup([atom for atom in Res_dn.atoms if atom.name == resT2[options.res]])

############################################################################


####################### SEARCH LOOP ########################################

n = 10

GX = np.zeros((u.trajectory.n_frames, n+1))
GY = np.zeros((u.trajectory.n_frames, n+1))

DS = np.zeros((u.trajectory.n_frames, 2, len(Lf_up), 4))
PS = np.zeros((u.trajectory.n_frames, 2, len(Lf_up), 3))

HS    = np.zeros((u.trajectory.n_frames, 2, n, n))
DmatS = np.zeros((u.trajectory.n_frames, 2, n, n, 3, 3)) 

def Compute(frame, Leaflet_up, Leaflet_down, TA_up, TA_dn, TB_up, TB_dn, n):

    Leaflet_up.universe.trajectory[frame]
    Leaflet_down.universe.trajectory[frame]

    TA_up.universe.trajectory[frame]
    TA_dn.universe.trajectory[frame]

    TB_up.universe.trajectory[frame]
    TB_dn.universe.trajectory[frame]

    Pup = Leaflet_up.positions
    Pdn = Leaflet_down.positions
   
    maxX = np.max([np.max(Pup[:,0]), np.max(Pdn[:,0])])
    maxY = np.max([np.max(Pup[:,1]), np.max(Pdn[:,1])])
    
    gX, gY = np.linspace(0, maxX, n+1), np.linspace(0, maxY, n+1)
    GX[frame] = gX
    GY[frame] = gY
    
    DT_up = Director(Lf_up, TA_up, TB_up)
    DT_dn = Director(Lf_dn, TA_dn, TB_dn)
    
    r_up = rot.from_matrix(DT_up)
    r_dn = rot.from_matrix(DT_dn)
    q_up = r_up.as_quat()
    q_dn = r_dn.as_quat()

    H_up, D_up = Search(n, Lf_up.positions, DT_up, gX, gY)
    H_dn, D_dn = Search(n, Lf_dn.positions, DT_dn, gX, gY)
    
    ##############################
    
    '''
    We will ensure that topography arrays dont have any empty values as this
    would break following computations. 
    
    This is achieved by substituting the 0. values with a linear interpolation
    of the contiguous bins.
    
    '''
    
    Interpolate(H_up, n)
    Interpolate(H_dn, n)
    
    #############################

    return gX, gY, Pup, Pdn, q_up, q_dn, H_up, H_dn, D_up, D_dn

#############################################################################

if __name__=='__main__':

    import multiprocessing
    from functools import partial

    run_frames = partial(Compute, Leaflet_up = Lf_up, Leaflet_down = Lf_dn, TA_up = TA_up, TA_dn = TA_dn, TB_up = TB_up, TB_dn = TB_dn, n = n)

    Frames = np.arange(u.trajectory.n_frames)

    pool = multiprocessing.Pool(options.prc)

    Result = pool.map(run_frames, Frames)

    GX += np.array([frame[0] for frame in Result])
    GY += np.array([frame[1] for frame in Result])

    PS[:,0] += np.array([frame[2] for frame in Result])
    PS[:,1] += np.array([frame[3] for frame in Result])

    DS[:,0] += np.array([frame[4] for frame in Result])
    DS[:,1] += np.array([frame[5] for frame in Result])

    HS[:,0] += np.array([frame[6] for frame in Result])
    HS[:,1] += np.array([frame[7] for frame in Result])

    DmatS[:,0] += np.array([frame[8] for frame in Result])
    DmatS[:,1] += np.array([frame[9] for frame in Result])

################### OUTPUT ##################################################

with open('gridX.pickle', 'wb') as f:
    pickle.dump(GX,f, protocol = 4)
f.close()

with open('gridY.pickle', 'wb') as f:
    pickle.dump(GY,f, protocol = 4)
f.close()

with open('topog.pickle', 'wb') as f:
    pickle.dump(HS,f, protocol = 4)
f.close()

with open('direc.pickle', 'wb') as f:
    pickle.dump(DS,f, protocol = 4)
f.close()

with open('posit.pickle', 'wb') as f:
    pickle.dump(PS,f, protocol = 4)
f.close()

with open('d_mat.pickle', 'wb') as f:
    pickle.dump(DmatS,f, protocol = 4)
f.close()
