#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 14:04:08 2022

@author: alexandre
"""

import numpy as np
import pickle
from optparse import OptionParser
from math import cos

parser = OptionParser(usage       = 'python thickness.py -k Pickles -o output -n <n_proc>',
                      prog        = 'Thickness',
                      description = 'This program calculates the thickness of'
                                    ' the membrane taking into account its curvature.')

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
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
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

def GetThickness(h_up, h_dn, lrs_up, lrs_dn):
    
    '''
    Calculate the thickness of a membrane taking into account its curvature.

    Parameters
    ----------
    h_up : NumPy array
        Array (X,Y) containing the surface of the upper leaflet.
    h_dn : NumPy array
        Array (X,Y) containing the surface of the lower leaflet.
    lrs_up : NumPy array
        Array (X,Y,3,3) containing the LRS on the upper leaflet.
    lrs_dn : NumPy array
        Array (X,Y,3,3) containing the LRS on the lower leaflet.

    Returns
    -------
    thick : NumPy array
        Array (X,Y) containing the thickness of the membrane.

    '''

    thick = np.zeros_like(h_up)
    Zup = np.array([0,0,1])
    Zdn = np.array([0,0,-1])

    for i in range(h_up.shape[0]):
        for j in range(h_up.shape[1]):
            
            zup = h_up[i,j]
            zdn = h_dn[i,j]
            
            nup = lrs_up[i,j,2]
            ndn = -1*lrs_dn[i,j,2]

            aup = angle_between(nup, Zup)
            adn = angle_between(ndn, Zdn)
            
            p2up = zup/cos(aup) # Change of basis to LRS
            p2dn = zdn/cos(adn)
            
            thick[i,j] += p2up - p2dn     # Thickness is Z'_up - Z'_down
    
    return thick

################# LOAD AUXILIARY TOPOGRAPHY AND LRS ARRAYS ##################

with open('./%s/topog.pickle' % options.kle, 'rb') as f:
    H = pickle.load(f)
f.close()

with open('./%s/lrss.pickle' % options.kle, 'rb') as f:
    LRS = pickle.load(f)
f.close()

with open('./%s/gridX.pickle' % options.kle, 'rb') as f:
    gX = pickle.load(f)
f.close()

with open('./%s/gridY.pickle' % options.kle, 'rb') as f:
    gY = pickle.load(f)
f.close()

#############################################################################


################# SEARCH LOOP ###############################################

THI = np.zeros((H.shape[0], H.shape[2], H.shape[3]))

def Compute( frame ):
    #print('Analyzed %i %% Frames' % (100*(frame+1)/l), end = '\r')
    
    gx = gX[frame]
    gy = gY[frame]
    
    x = 0.5*(gx[1:]+gx[:-1])
    y = 0.5*(gy[1:]+gy[:-1])
    
    h_up = H[frame,0]
    h_dn = H[frame,1]
    
    lrs_up = LRS[frame,0]
    lrs_dn = LRS[frame,1]
    
    thick = GetThickness(h_up, h_dn, lrs_up, lrs_dn)

    return  thick 

if __name__ == '__main__':
    
    import multiprocessing as mp
    pool = mp.Pool(options.prc)
    Result = pool.map( Compute, range( H.shape[0] ) )
    THI += np.array( Result )

#############################################################################


################# OUTPUT ####################################################

with open('%s.pickle' % options.out, 'wb') as f:
    pickle.dump(THI,f, protocol = 4)
f.close()
