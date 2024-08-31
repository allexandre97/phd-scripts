#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 14:33:40 2022

@author: alexandre
"""

import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder as lf
import numpy as np
import pyvista as pv
from optparse import OptionParser
import time

parser = OptionParser(usage       = 'python area_per_lipid.py -p input.pdb -x input.xtc -r <resolution> -o output -n <n_proc>',
                      prog        = 'AreaPerLipid',
                      description = 'This program calculates the area per lipid'
                                    ' from simulation data performing a triangulation'
                                    ' over the surface of the membrane.')

parser.add_option('-p', '--pdb',
                  action='store', type = 'string',
                  help = 'Path to input structure file')

parser.add_option('-x', '--xtc',
                  action='store', type = 'string',
                  help = 'Path to input trajectory file')

parser.add_option('-r', '--res',
                  action='store', type = 'string',
                  help = 'Resolution of the simulation (AA or CG)')

parser.add_option('-o', '--out',
                  action='store', type = 'string',
                  help = 'Output pickled file')

parser.add_option('-n', '--prc',
                  action='store', type = 'int',
                  help = 'Number of Processors, defaults to 6',
                  default = 6)

options, arguments = parser.parse_args()

def AreaPerFace(Faces, Points):
    '''
    Compute the area of each of the faces of the Triangulated Surface

    Parameters
    ----------
    Faces : NumPy array
        Padded array containing the indices of the vertices defining each face.
    Points : NumPy array
        Array containing the positions of each vertex defining each face.

    Returns
    -------
    Area : NumPy array
        Area of each face.

    '''
    
    i = 0 # Auxiliary ints defining the start and end of each face
    f = 4 # in the index array. Padding is constant if surface is a Delaunay
          # triangulation. THIS FUNCTION WON'T WORK OTHERWISE!!!!
          
    Area = []
    while f <= Faces.shape[0]:
        face = Faces[i+1:f]   # Get the indices defining the face
        points = Points[face] # Get corresponding points 
        
        ''' 
        Area of a 3D triangle (ABC) can be caclulated as follows:
        First define two of the vertices, i.e. AC and AB.
        Area is 0.5*|AB x AC|
        '''
        
        v1 = points[0] - points[1]
        v2 = points[0] - points[2]
        cross = np.cross(v1,v2)
        Cross = np.sqrt(cross[0]**2+cross[1]**2+cross[2]**2)
        area = 0.5*Cross
        Area.append(area)
        
        i = f
        f = i+4
    
    Area = np.array(Area)
    return Area

def VertexArea(Faces, faces, points, area, p):
    '''
    Compute the surrounding area of a vertex as the sum of the corner areas 
    incident upon the vertex.

    Parameters
    ----------
    Faces : NumPy array
        Array containing the indices of the vertex composing each face of the 
        triangulated surface.
    faces : NumPy array
        Array containing the indices of the faces surrounding the vertex of
        interest.
    points : NumPy array
        Array containing the 3D positions of the points defining the surface.
    area : NumPy array
        Array containing the area of each face.
    p : Int
        index of the vertex of interest.

    Returns
    -------
    Suma : Float
        Area around the vertex of interest.

    '''
    
    Suma = 0
    for f in faces:

        Af = area[f] # Get the area of the face
        
        ids = Faces[f*4+1:f*4+4] # Get the vertex indices defining the face
        
        k = np.where(ids==p)[0][0] # Get the index of the vertex of interest
        
        k1 = (k+1) % 3 # Mod operator applied to ids of the vertices 
        k2 = (k+2) % 3 # as to make k the "first" index
        
        '''
        The computation of the area follows the procedure described in:

        MemSurfer: A Tool for Robust Computation and Characterization of Curved Membranes, 
        https://doi.org/10.1021/acs.jctc.9b00453
        
        As a quick explanation, the corner areas of a triangle can be compute using the
        barycentric coordinates of the circumcenter of the triangle. It follows:

        a point P in a triangle (ABC) can be expressed by means of barycentric coordinates such as:

        P = mu_A * A + mu_B * B + mu_C * C; 
        mu_A + mu_B + mu_C = 1. --> These are the barycentric coordinates

        When P is the circumcenter, the barycentric coordinates are computed using the length
        of the sides of the triangle. Let:

        e_0 = C - B
        e_1 = A - C --> These are the sides of the triangle
        e_2 = B - A

        l_k = | e_k | --> The length of the kth side

        Then,

        mu_k = l_k * (l_(k+1 mod 3) + l_(k+2 mod 3) - l_k) --> kth barycentric coordinate.

        Note the mod operator, prevents overflow of index taking advantage of
        closed shape of triangle (Vertex 4 = Vertex 1)

        And finally, the corner area is given by

        CornerArea = AreaTriangle * (mu_(k+1 mod 3) + mu_(k+2 mod 3))/2 

        Where k is the index of the vertex of interest.
        '''
        
        e0 = points[ids[2]] - points[ids[1]]
        e1 = points[ids[0]] - points[ids[2]] # Compute the side vectors of the face
        e2 = points[ids[1]] - points[ids[0]]
        
        l0 = np.sqrt(e0[0]**2+e0[1]**2+e0[2]**2)
        l1 = np.sqrt(e1[0]**2+e1[1]**2+e1[2]**2) # Compute the length of the sides
        l2 = np.sqrt(e2[0]**2+e2[1]**2+e2[2]**2)
        
        ls = [l0,l1,l2]
        
        mu0 = ls[0]*(ls[1] + ls[2] - ls[0])
        mu1 = ls[1]*(ls[2] + ls[0] - ls[1]) # Compute the barycentric coordinates of circumcenter
        mu2 = ls[2]*(ls[0] + ls[1] - ls[2])
        
        mus  = np.array([mu0, mu1, mu2]) # Normalize mu
        mus /= (mus[0]+mus[1]+mus[2])
        
        Suma += Af*((mus[k1] + mus[k2])/2) # Add to total area
        
    return Suma

def AreaPerVertex(Points, Faces, area):
    '''
    Iterates over the vertices, finds the adjacent faces and calculates the
    area per vertex based on a barycentric criterion.

    Parameters
    ----------
    Points : NumPy array
        Array containing the 3D positions of the points defining the surface.
    Faces : NumPy array
        Array containing the indices of the vertex composing each face of the 
        triangulated surface.
    area : NumPy array
        Array containing the area of each face.

    Returns
    -------
    apl : NumPy array
        array containing the area corresponding to each vertex (lipid).

    '''
    
    apl = np.zeros((Points.shape[0])) # Initialize output array
    

    for p in range(Points.shape[0]):
        id_faces = np.where(Faces==p)[0] # Find where the point appears
                                         # in the faces array
        faces = []
        '''
        These conditionals take into account that the faces array is padded.
        Because faces are always triangles, padding is always 3, so when we
        look at point with idx = 3, we have to be carful to not confuse 
        the actual appearance of the point in a face with the padding.
        '''
        if p != 3:
            for i in id_faces:
                faces.append(int(np.floor(i/4))) # Get an index of face. Array
                                                 # is 3-padded, so each face is 
                                                 # described by 1+3 numbers: 
                                                 # padding + vertices
            faces = np.array(faces)
            VArea  = VertexArea(Faces, faces, Points, area, p)
            apl[p] = VArea

        else:
            for n, i in enumerate(id_faces[:-1]):
                if abs(i - id_faces[n+1]) != 4: # Ensure that the 3 we find actually
                                                # represents a point

                    faces.append(int(np.floor(i/4)))
            faces = np.array(faces)
            faces = np.unique(faces)
            VArea  = VertexArea(Faces, faces, Points, area, p)
            apl[p] = VArea
            
    return apl

################ LOAD GROMACS TRAJECTORY AND DIVIDE LEAFLETS ################

u = mda.Universe(options.pdb, options.xtc, in_memory = False)

res = {'AA':'P', 'CG':'PO4'}

P = u.select_atoms('name %s' % res[options.res], periodic = True, updating = True)

leafs = lf(u, P, pbc = True)

Lf_up = leafs.groups(0)
Lf_dn = leafs.groups(1)

#############################################################################


#################### SEARCH LOOP ############################################

APL = np.zeros((u.trajectory.n_frames, 2, len(Lf_up)))

def Compute(frame, Lf_up, Lf_dn):

    #print(frame)
    
    Lf_up.universe.trajectory[frame]
    Lf_dn.universe.trajectory[frame]

    Pup = Lf_up.positions
    Pdn = Lf_dn.positions
    
    cloudup = pv.PolyData(Pup)
    clouddn = pv.PolyData(Pdn)
    
    surfup = cloudup.delaunay_2d()
    surfdn = clouddn.delaunay_2d()
    
    faces_up = surfup.faces
    faces_dn = surfdn.faces
    
    points_up = surfup.points
    points_dn = surfdn.points
    
    Area_up = AreaPerFace(faces_up, points_up)
    Area_dn = AreaPerFace(faces_dn, points_dn)
    
    apl_up = AreaPerVertex(points_up, faces_up, Area_up)
    apl_dn = AreaPerVertex(points_dn, faces_dn, Area_dn)

    #time.sleep(0.005)

    return apl_up, apl_dn

#############################################################################

if __name__=='__main__':

    import multiprocessing
    from functools import partial

    run_frames = partial(Compute, Lf_up = Lf_up, Lf_dn = Lf_dn)

    Frames = np.arange(u.trajectory.n_frames)

    pool = multiprocessing.Pool(options.prc)

    print('Computing')

    Result = pool.map(run_frames, Frames)

    print('Done!')

    APL[:,0] = np.array([frame[0] for frame in Result])
    APL[:,1] = np.array([frame[1] for frame in Result])

################## OUTPUT ###################################################

import pickle

with open('%s.pickle' % options.out, 'wb') as f:
    pickle.dump(APL,f, protocol = 4)
f.close()
