#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 16:22:57 2022

@author: alexandre
"""

import pickle
import numpy as np
from scipy.spatial.transform import Rotation as rot
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter as sf

def forceAspect(ax,aspect):
    '''
    Forces square aspect in a matplotlib figure.
    
    Parameters
    ----------
    ax : Matplotlib axis object.
        The subplot to apply aspect to.
    aspect : float
        Aspect ratio of subplot, 1 is square.

    Returns
    -------
    None.

    '''
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def bins(X):
    '''
    Generate bins for a histogram following the Friedman-Diaconis rule.

    Parameters
    ----------
    X : numpy array
        The variable to calculate the histogram for.

    Returns
    -------
    Bins : numpy array
        Binning for histogram.

    '''
    q25, q75 = np.percentile(X, [25, 75])
    iqr = q75 - q25
    binwidth = (2*iqr)/(len(X)**(1/3))*(1-0.2)
    Bins = np.arange(min(X), max(X) + binwidth, binwidth)
    return Bins

def unit_vector(vector):
    '''
    Normalizes a vector so it is of unit length.

    Parameters
    ----------
    vector : numpy array
        N-dimensional vector.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    return vector / np.linalg.norm(vector)

def cosine_between(v1, v2):
    '''
    Calculate the projection of vector 2 over vector 1.

    Parameters
    ----------
    v1 : numpy array
        First N-dimensional vector.
    v2 : numpy array
        Second N-dimensional vector.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)

def DirectorTilt(gx, gy, position, lrs, director):
    '''
    Looks for the corresponding LRS of each lipid director and calculates its
    tilt as the cosine of the angle between the director and the normal.

    Parameters
    ----------
    gx : numpy array
        Grid over the membrane over the X dimension.
    gy : numpy array
        Grid over the membrane over the Y dimension.
    position : numpy array
        Array with the positions of the phosphate groups, used to check the 
        correspondent LRS.
    lrs : numpy array
        Local Reference System in each cell grid.
    director : numpy array
        Lipid director of rach lipid..

    Returns
    -------
    Cos : numpy array
        Projections of the lipid director over the normals.

    '''
    Cos = np.zeros((position.shape[0]))
    for i, p in enumerate(position):
        setx = gx < p[0] # Gridx smaller than position X
        sety = gy < p[1] # Gridy smaller than position Y
        idx = np.where(setx==False)[0][0]-1 # Get indices of grid where it
        idy = np.where(sety==False)[0][0]-1 # is no longer smaller than pX & pY
        n = lrs[idx,idy,2] # Normal vector at selected grid point
        v = director[i]  # Lipid director vector
        cos = cosine_between(n, v) # Calculate cosine between vectors
        Cos[i] += cos # Store angle
    return Cos

with open('properties.xvg', 'r') as f:
    lines = f.readlines()
    t = 4*np.array([float(_.split()[0]) for _ in lines[45:]])/1000
    Pxx = np.array([float(_.split()[13]) for _ in lines[45:]])
    Pxy = np.array([float(_.split()[14]) for _ in lines[45:]])
    Pxz = np.array([float(_.split()[15]) for _ in lines[45:]])
    Pyx = np.array([float(_.split()[16]) for _ in lines[45:]])
    Pyy = np.array([float(_.split()[17]) for _ in lines[45:]])
    Pyz = np.array([float(_.split()[18]) for _ in lines[45:]])
    Pzx = np.array([float(_.split()[19]) for _ in lines[45:]])
    Pzy = np.array([float(_.split()[20]) for _ in lines[45:]])
    Pzz = np.array([float(_.split()[21]) for _ in lines[45:]])
    T = np.array([float(_.split()[22]) for _ in lines[45:]])
f.close()

PRE = np.zeros((t.shape[0],3,3))
PRE[:,0,0] += Pxx
PRE[:,0,1] += Pxy
PRE[:,0,2] += Pxz
PRE[:,1,0] += Pyx
PRE[:,1,1] += Pyy
PRE[:,1,2] += Pyz
PRE[:,2,0] += Pzx
PRE[:,2,1] += Pzy
PRE[:,2,2] += Pzz

PR = sf(PRE, 80*50+1, 2, axis = 0)

pmin = np.min(PR[:,2,2])*(1 - 0.1)
pmax = np.max(PR[:,2,2])*(1 + 0.1)

with open('direc.pickle', 'rb') as f:
    DIR = pickle.load(f)
f.close()

with open('posit.pickle', 'rb') as f:
    POS = pickle.load(f)
f.close()

with open('gridX.pickle', 'rb') as f:
    GX = pickle.load(f)
f.close()

with open('gridY.pickle', 'rb') as f:
    GY = pickle.load(f)
f.close()

with open('lrss.pickle', 'rb') as f:
    LRS = pickle.load(f)
f.close()


dt = 20

w = int(np.ceil(5000/dt))

Hist = []
Mean = []

bns = np.linspace(-1, -0.5, 500)
    
def ComputeTilt(frame, POS, LRS, DIR):
    
    gx = GX[frame]
    gy = GY[frame]

    position_up = POS[frame,0]
    position_dn = POS[frame,1]
    
    director_up = rot.from_quat(DIR[frame,0]).as_matrix()[:,0]
    director_dn = rot.from_quat(DIR[frame,1]).as_matrix()[:,0]
    
    lrs_up = LRS[frame,0]
    lrs_dn = LRS[frame,1]
    
    Cos_up = DirectorTilt(gx, gy, position_up, lrs_up, director_up)
    Cos_dn = DirectorTilt(gx, gy, position_dn, lrs_dn, director_dn)
    
    return [Cos_up, Cos_dn]

def ComputeHistogram(frame, window, COS):
    
    cos_up   = COS[frame:frame+window,0].flatten()    
    m_up     = np.mean(cos_up)
    s_up     = np.std(cos_up)
    h_up, _  = np.histogram(cos_up, bins = bns, density = True)

    cos_dn   = -1*COS[frame:frame+window,1].flatten()
    m_dn     = np.mean(cos_dn)
    s_dn     = np.std(cos_dn)
    h_dn, _  = np.histogram(cos_dn, bins = bns, density = True)
    
    return [h_up, h_dn], [[m_up, s_up], [m_dn, s_dn]]

if __name__ == '__main__':
    
    import multiprocessing as mp
    from functools import partial
    
    run_tilt = partial(ComputeTilt,
                       POS = POS,
                       LRS = LRS,
                       DIR = DIR)
    
    pool = mp.Pool(10)
    
    Frames = np.arange(POS.shape[0])
    
    print('Calculating Tilt')
    
    COS = pool.map(run_tilt, Frames)
    
    COS = np.array(COS)
    
    run_hist = partial(ComputeHistogram,
                       window = w,
                       COS = COS)
    
    Frames = np.arange(COS.shape[0] - w)
    
    print('Calculating Histogram')
    
    Result = pool.map(run_hist, Frames)
    
    Hist = np.array([frame[0] for frame in Result])
    Mean = np.array([frame[1] for frame in Result])
    

with open('tails_dist.pickle', 'wb') as f:
    pickle.dump(Hist, f)
f.close()

with open('tails_mean.pickle', 'wb') as f:
    pickle.dump(Mean, f)
f.close()

T = 0.02*np.linspace(0, Hist.shape[0], Hist.shape[0])

a = (t[-1]-t[0])/90

Min, Max = np.min([Hist[:,0].min(), Hist[:,1].min()]), np.max([Hist[:,0].max(), Hist[:,1].max()])*(1-0.75)

fig, ax = plt.subplots(1, 2, figsize = (12,5))

ax[0].set_title('Upper Leaflet', fontsize = 22, fontweight = 'bold', pad = 60)
im1 = ax[0].imshow(np.log(Hist[:,0].T + 1e-9),
                   cmap = 'inferno',
                   origin = 'lower',
                   extent = [t[0], t[-1], bns[0], bns[-1]],
                   vmin = Min, vmax = Max,
                   rasterized = True)
ax[0].plot(T, Mean[:,0,0], color = 'k')
ax[0].set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold')
ax[0].set_ylabel('$\mathbf{\Theta (º)}$', fontsize = 20, fontweight = 'bold')
ax[0].tick_params(labelsize=18)
ax[0].set_ylim(bns[0]-0.1, bns[-1])
forceAspect(ax[0], aspect=1.0)
ax2 = ax[0].twinx()
ax2.plot(t, PR[:,2,2], color = 'w')
ax2.set_ylabel('Pr. (bar)', fontsize = 12, fontweight = 'bold', labelpad = -15)
ax2.set_ylim(pmin-10, pmax)
ticks = np.linspace(pmin-10, pmax, 4)
ax2.set_yticks(ticks)
ax2.set_yticklabels(["{:2.0f}".format(int(i)) for i in ticks], fontsize='14')

v1 = np.linspace(Min, Max, 4, endpoint=True)
cb = fig.colorbar(im1, ticks=v1, ax = ax2, shrink = 1, location = 'top')
cb.set_label(r'$\mathbf{P(\Theta)}$',
                fontsize = 11,
                fontweight = 'bold',
                rotation = 0,
                labelpad = -15)
cb.ax.set_xticklabels(["{:4.2f}".format(i) for i in v1], fontsize='14')

ax[1].set_title('Lower Leaflet', fontsize = 22, fontweight = 'bold', pad = 60)
im2 = ax[1].imshow(np.log(Hist[:,1].T + 1e-9),
                   cmap = 'inferno',
                   origin = 'lower',
                   extent = [t[0], t[-1], bns[0], bns[-1]],
                   vmin = Min, vmax = Max,
                   rasterized = True)
ax[1].plot(T, Mean[:,1,0], color = 'k')
ax[1].set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold')
ax[1].set_ylabel('$\mathbf{\Theta (º)}$', fontsize = 20, fontweight = 'bold')
ax[1].tick_params(labelsize=18)
ax[1].set_ylim(bns[0]-0.1, bns[-1])
forceAspect(ax[1], aspect=1.0)
ax2 = ax[1].twinx()
ax2.plot(t, PR[:,2,2], color = 'w')
ax2.set_ylabel('Pr. (bar)', fontsize = 12, fontweight = 'bold', labelpad = -15)
ax2.set_ylim(pmin-10, pmax)
ticks = np.linspace(pmin-10, pmax, 4)
ax2.set_yticks(ticks)
ax2.set_yticklabels(["{:2.0f}".format(int(i)) for i in ticks], fontsize='14')

v1 = np.linspace(Min, Max, 4, endpoint=True)
cb2 = fig.colorbar(im2, ticks=v1, ax = ax2, shrink = 1, location = 'top')
cb2.set_label(r'$\mathbf{P(\Theta)}$',
                fontsize = 11,
                fontweight = 'bold',
                rotation = 0,
                labelpad = -15)
cb2.ax.set_xticklabels(["{:4.2f}".format(i) for i in v1], fontsize='14')

fig.tight_layout()

plt.savefig('tailtilt.svg', dpi = 300)