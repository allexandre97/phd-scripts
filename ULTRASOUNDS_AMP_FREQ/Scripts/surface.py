#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 12:45:38 2022

@author: alexandre
"""

import numpy as np
import pickle
from scipy.interpolate import RectBivariateSpline as rbs
from optparse import OptionParser

parser = OptionParser(usage       = 'python surface.py -k <path_to_pickle_files> -n <n_proc> ',
                      prog        = 'Surface',
                      description = 'This program calculates several important.'
                                    ' structural and geometrical parameters of'
                                    ' the surface of the memebrane, such as its'
                                    ' curvature, the Local Reference System and'
                                    ' the normals.')

parser.add_option('-k', '--kle',
                  action='store', type = 'string',
                  help = 'Path to Pickles folder')

parser.add_option('-n', '--prc',
                  action='store', type = 'int',
                  help = 'Number of Processors, defaults to 6',
                  default = 6)

options, arguments = parser.parse_args()


def Gradient(H, dx):
    '''
    Calculate the gradient of a 2D scalar field. This is calculated via a 
    centered finite element method, and accounting for Periodic Boundary 
    Conditions.

    Parameters
    ----------
    H : NumPy array
        2D array of the scalar field for which the gradient will be calculated.
    dx : Float
        Differential increment, i.e grid spacing.

    Returns
    -------
    G : NumPy array
        Array (X,Y,2) containing the gradient.

    '''
    dx2 = 2*dx

    G = np.zeros((H.shape[0],H.shape[1],2))
    G[0,:,0]       = (H[1,:]-H[-1,:])/dx2
    G[0,1:-1,1]    = (H[0,2:]-H[0,:-2])/dx2
    G[-1,:,0]      = (H[0,:]-H[-2,:])/dx2
    G[-1,1:-1,1]   = (H[-1,2:]-H[-1,:-2])/dx2
    G[1:-1,0,0]    = (H[2:,0]-H[:-2,0])/dx2
    G[:,0,1]       = (H[:,1]-H[:,-1])/dx2
    G[1:-1,-1,0]   = (H[2:,-1]-H[:-2,-1])/dx2
    G[:,-1,1]      = (H[:,0]-H[:,-2])/dx2
    G[1:-1,1:-1,0] = (H[2:,1:-1]-H[:-2,1:-1])/dx2
    G[1:-1,1:-1,1] = (H[1:-1,2:]-H[1:-1,:-2])/dx2
    return G

def Jacobian(f, dx):
    '''
    Compute the Jacobian (2D) of a vector valued function. This is done through
    the gradient of each one of the components of the vector funcion, so it 
    inherits the finite element method of the Gradient function.

    Parameters
    ----------
    f : NumPy array
        Input vector field.
    dx : Float
        Differential increment, i.e grid spacing.

    Returns
    -------
    Jacob : NumPy array
        Array (X,Y,2,2) containig the Jacobian at each point of the field.

    '''
    Jacob = np.zeros((f.shape[0], f.shape[1], 2, 2))
    for c in range(2):
        # Note that jacobian matrix is the matrix of gradients
        # of a vector valued function
        Jacob[:,:,c] += Gradient(f[:,:,c], dx)
    return Jacob

def LTV(G, g):
    '''
    Calculate the Local Tangent Vectors (plus normals!) of a surface.

    Parameters
    ----------
    G : NumPy array
        Array (X,Y,2) containing the Gradient of the scalar field (i.e. elevation)
        defining the surface.
    g : NumPy array
        Array (X,Y) containing the metric tensor determinant at each point in
        the surface. Important for membranes that deviate a lot from planarity.

    Returns
    -------
    ltv : NumPy array
        Array (X,Y,3,3) containig the Local Tangent Vectors to the surface plus
        its normals at each point.

    '''
    
    ltv = np.zeros((G.shape[0], G.shape[1], 3, 3)) # Array for Local Tangent vectors
                                                   # Third vector is for normals
    
    '''
    The framework of this analysis is within Monge Gauge theory, where a surface 
    is parametrized as the (X, Y) points and an elevation function (h) over or
    below the XY plane.
    
    Thus, the Local Tangent vectors can be defined as e1 = (1, 0, dh/dx) and 
    e2 = (0, 1, dh/dy).
    
    Then, the normals can be computed as  (e1 x e2)/sqrt(g) where g is the 
    determinant of the metric tensor.
    
    For a very in good review on why this is true and how it works check --->
    Fluid lipid membranes: From differential geometry to curvature stresses 
    by Markus Deserno
    '''
    
    ltv[:,:,0,0] += 1
    ltv[:,:,0,2] += G[:,:,0]
    ltv[:,:,1,1] += 1
    ltv[:,:,1,2] += G[:,:,1]
    for i in range(ltv.shape[0]):
        for j in range(ltv.shape[1]):
            ltv[i,j,2] += np.cross(ltv[i,j,0], ltv[i,j,1])/np.sqrt(g[i,j])
    return ltv

def Curv(G, g, dx):
    '''
    Compute the total curvature of a surface. Done through the Jacobian of the
    Gradient (the Hessian!) of the elevation defining function
    parametrizing the surface. 

    Parameters
    ----------
    G : NumPy array
        Array (X,Y,2) containing the Gradient of the scalar field (i.e. elevation)
        defining the surface.
    g : NumPy array
        Array (X,Y) containing the determinant of the metric tensor at each point in
        the surface. Important for membranes that deviate a lot from planarity.
    dx : Float
        Differential increment, i.e grid spacing.
        
    Returns
    -------
    K : NumPy array
        Array (X,Y) containing the total curvature of the surface at each point.

    '''
    
    '''
    Again, within Monge Gauge we can compute the total curvature of a surface
    as -Tr{Jac[G(h)/sqrt(g)]} = -Tr[Hess(h)] where:
        
    ---> Tr is the trace of a matrix
    ---> Jac is the Jacobian of a vector field
    ---> h is the elevation function 
    ---> g is the determinant of the metric tensor
    ---> Hess is the Hessian of a scalar field.
    '''
    
    tmp = np.zeros_like(G)
    for i in range(tmp.shape[0]):
        for j in range(tmp.shape[1]):
            tmp[i,j] += G[i,j]/np.sqrt(g[i,j])  # First step is to compute the quotient
                                                # of the gradient of the elevation
                                                # and the metric determinant
                                                
    tmp = Jacobian(tmp, dx) # Compute the Jacobian on the previous quotient
    K = np.zeros((G.shape[0], G.shape[1]))
    for i in range(K.shape[0]):
        for j in range(K.shape[1]):
            K[i,j] = -np.trace(tmp[i,j]) # Curvature is the trace 
                                         # of the previous Jacobian
    return K

def LRS(K, GH, LTV, dx):
    '''
    Compute the Local Reference Systems over a surface. LRS is defined as 
    the Isocurvature, Anisocurvature an Normal vector fields.
    
    Parameters
    ----------
    K : NumPy array
        Array (X,Y) containing the total curvature of the system at each point.
    GH : NumPy array
        Array (X,Y,2) containing the Gradient of the elevation function of the 
        surface at each point.
    LTV : NumPy array
        Array (X,Y,3,3) containing the Local Tangent Vectors (plus normals) of
        the surface at each point.
    dx : Float
        Differential increment, i.e grid spacing.

    Returns
    -------
    LRS : NumPy array
        Array (X,Y,3,3) containing the Local Reference System at each point.

    '''
    
    '''
    By definition, LRS will be the directions of minimum (iso) curvature,
    maximum (aniso) curvature and the normals on the surface. 
    
    Thus the first step is to calculate the Gradient of the total curvature (K)
    scalar field. However, the gradient will only tell about the X and Y
    directions. To get the Z component of this field, we can use the gradient 
    of the topography (eleveation), as the gradient of the curvature must 
    follow the surface. 
    
    Hence, the Z component of the curvature (K) gradient can be calculated as the 
    dot product of said Gradient in XY and the Gradient of the topography (h).
    
    ---> GK = (dK/dx, dK/dy, [dK/dx,dK/dy]·[dh/dx, dh/dy])
    
    However, there is no guarantee that this vector field will be orthogonal to 
    the normals (N). To solve this, we define the Anisocurvature as
    
    ---> Anis. = (N x GK) / |N x GK|, which ensures orthonormality.
    
    Now, we can define the Isocurvature as:
    
    ---> Isoc. = (N x Anis.) / |N x Anis|
    
    '''
    
    GK = Gradient(K,dx) #First take the gradient of the curvature
    lrs = np.zeros((K.shape[0], K.shape[1], 3, 3))
    for i in range(lrs.shape[0]):
        for j in range(lrs.shape[1]):
            z = GK[i,j].dot(GH[i,j])  # This gives the Z component of the Curvature gradient
            v = np.array([GK[i,j,0],GK[i,j,1],z])
            n = LTV[i,j,2] # Normal Vector
            tmp = np.cross(n, v)
            tmp /= np.sqrt(tmp[0]**2+tmp[1]**2+tmp[2]**2)
            lrs[i,j,0] += tmp
            tmp = np.cross(n,tmp)
            tmp /= np.sqrt(tmp[0]**2+tmp[1]**2+tmp[2]**2)
            lrs[i,j,1] += tmp
            lrs[i,j,2] += n
    
    return lrs

################## LOAD AUXILIARY TOPOGRAPHY AND GRID ARRAYS ################

with open('./%s/gridX.pickle' % options.kle, 'rb') as f:
    gX = pickle.load(f)
f.close()

with open('./%s/gridY.pickle' % options.kle, 'rb') as f:
    gY = pickle.load(f)
f.close()

with open('./%s/topog.pickle' % options.kle, 'rb') as f:
    H = pickle.load(f)
f.close()

#############################################################################


################ SEARCH LOOP ################################################

LTVS = np.zeros((H.shape[0], H.shape[1], H.shape[2], H.shape[3], 3, 3))
LRSS = np.zeros_like(LTVS)
KS   = np.zeros((H.shape[0], H.shape[1], H.shape[2], H.shape[3]))

def Compute(frame, gX, gY, H):
     
    gx, gy = gX[frame], gY[frame]
    X = 0.5*(gx[1:]+gx[:-1])
    Y = 0.5*(gy[1:]+gy[:-1])
    
    ################
    h_up = H[frame,0]

    h_dn = H[frame,1]
    #####################
     
    dx = X[1]-X[0]
    
    Gh_up = Gradient(h_up, dx)
    Gh_dn = Gradient(h_dn, dx)
    
    g_up = np.ones((Gh_up.shape[0], Gh_up.shape[1])) #Array for metric determinant
    g_up += Gh_up[:,:,0]**2 + Gh_up[:,:,1]**2

    g_dn = np.ones((Gh_dn.shape[0], Gh_dn.shape[1])) #Array for metric determinant
    g_dn += Gh_dn[:,:,0]**2 + Gh_dn[:,:,1]**2
    
    LTV_up = LTV(Gh_up, g_up)  # Local Tangent Vectors
    LTV_dn = LTV(Gh_dn, g_up)
    
    K_up = Curv(Gh_up, g_up, dx) # Total curvatures
    K_dn = Curv(Gh_dn, g_dn, dx)
    
    LRS_up = LRS(K_up, Gh_up, LTV_up, dx) # Local Reference System
    LRS_dn = LRS(K_dn, Gh_dn, LTV_dn, dx)
    
    return LTV_up, LTV_dn, LRS_up, LRS_dn, K_up, K_dn
   
#############################################################################

if __name__=='__main__':

    import multiprocessing
    from functools import partial

    run_frames = partial(Compute, gX = gX, gY = gY, H = H)

    Frames = np.arange(H.shape[0])

    pool = multiprocessing.Pool(options.prc)

    Result = pool.map(run_frames, Frames)

    LTVS[:,0] += np.array([frame[0] for frame in Result])
    LTVS[:,1] += np.array([frame[1] for frame in Result])
    
    LRSS[:,0] += np.array([frame[2] for frame in Result])
    LRSS[:,1] += np.array([frame[3] for frame in Result])

    KS[:,0] += np.array([frame[4] for frame in Result])
    KS[:,1] += np.array([frame[5] for frame in Result])


######################### OUTPUT ############################################

with open('curv.pickle', 'wb') as f:
    pickle.dump(KS,f, protocol = 4)
f.close()

with open('ltvs.pickle', 'wb') as f:
    pickle.dump(LTVS, f, protocol = 4)
f.close()

with open('lrss.pickle','wb') as f:
    pickle.dump(LRSS,f, protocol = 4)
f.close()
