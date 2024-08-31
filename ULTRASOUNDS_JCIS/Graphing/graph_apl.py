#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 16:36:44 2022

@author: alexandre
"""

import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter as sf

def forceAspect(ax,aspect):
    '''Force the aspect of the graph to a desired value (1 is a square),
    useful for non-square images'''
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)


'''Load pressure tensor from properties xvg file'''
with open('pressure.xvg', 'r') as f:
    TP = np.loadtxt(f, comments = ('#', '@'))
f.close()

'''Store pressure tensor in one array'''
PRE = TP[:,1]
t   = 4*TP[:,0]

dt1 = t[1]-t[0]

'''Savitzsky-Golay filter on pressure to substract noise'''
PR = sf(PRE, 80*50+1, 2, axis = 0)

'''Minimum and maximum pressure values (plus a little bit of room around them)'''
pmin = np.min(PR)*(1+0.1)
pmax = np.max(PR)*(1+0.1)

'''Load Area Per Lipid pickle file'''
with open('apl.pickle', 'rb') as f:
    APL = pickle.load(f)
f.close()

bins = np.linspace(0, 1.5, 100)

'''Timestep parameters to define averaging window (w)
and Time array for mean plotting'''
dt = 20
w  = int(np.ceil(5000/dt))

dt2   = t[-1]/APL.shape[0]
idxi = int(1000*float(sys.argv[2])/dt2)
idxf = int(1000*float(sys.argv[3])/dt2)

APL = APL[idxi:idxf + w]

idxi = int(1000*float(sys.argv[2])/dt1)
idxf = int(1000*float(sys.argv[3])/dt1)

t = t[idxi:idxf]/1000
PR = PR[idxi:idxf]

def Compute(frame, window, APL):
    
    dummy = APL[frame:frame+w]
    
    allup = dummy[:,0].flatten()/100
    alldn = dummy[:,1].flatten()/100
    
    n_up, _ = np.histogram(allup, bins = bins, density = True)
    n_dn, _ = np.histogram(alldn, bins = bins, density = True)
    
    return [n_up, n_dn], [np.mean(allup), np.mean(alldn)]

if __name__ == '__main__':
    
    import os
    import multiprocessing as mp
    from functools import partial

    run = partial(Compute,
                  window = w,
                  APL = APL)

    Frames = np.arange(APL.shape[0]-w)

    pool = mp.Pool(int(sys.argv[1]))
    Result = pool.map(run, Frames)

    Image = np.array([frame[0] for frame in Result])
    Means = np.array([frame[1] for frame in Result])

import pickle

with open('apl_dist.pickle', 'wb') as f:
    pickle.dump(Image, f)
f.close()

with open('apl_mean.pickle', 'wb') as f:
    pickle.dump(Means, f)
f.close()

Time = np.linspace(t[0], t[-1], Means.shape[0])

'''Plotting!'''

Max, Min = np.max([np.max(Image[:,0]), np.max(Image[:,1])]), np.min([np.min(Image[:,0]), np.min(Image[:,0])])

fig, ax = plt.subplots(1,2,figsize=(12,5))

ax[0].set_title('Upper Leaflet', fontsize = 22, fontweight = 'bold', pad = 65)
im1 = ax[0].imshow(Image[:,0].T, origin = 'lower',
                   cmap = 'inferno',
                   extent=[t[0], t[-1], bins[0], bins[-1]],
                   vmin = Min, vmax = Max,
                   rasterized = True)
ax[0].plot(Time, Means[:,0], c = 'k')
ax[0].set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold')
ax[0].set_ylabel(r'$\mathbf{A.P.L. (\AA^{2})}$', fontsize = 20, fontweight = 'bold')
ax[0].tick_params(labelsize = 18)
forceAspect(ax[0], aspect=1.0)
ax2 = ax[0].twinx()
ax2.plot(t, PR, color = 'w')
ax2.set_ylabel('Pr. (bar)', fontsize = 12, fontweight = 'bold', labelpad = -15)
ax2.set_ylim(pmin, pmax)
ticks = np.linspace(pmin, pmax, 4)
ax2.set_yticks(ticks)
ax2.set_yticklabels(["{:2.0f}".format(int(i)) for i in ticks], fontsize='14')

v1 = np.linspace(Min, Max, 4, endpoint=True)
cb = fig.colorbar(im1, ticks=v1, ax = ax2, shrink = 1, location = 'top')
cb.set_label('P (A.P.L.)',
                fontsize = 11,
                fontweight = 'bold',
                rotation = 0,
                labelpad = 0)
cb.ax.set_xticklabels(["{:4.2f}".format(i) for i in v1], fontsize='14')

ax[1].set_title('Lower Leaflet', fontsize = 22, fontweight = 'bold', pad = 65)
im2 = ax[1].imshow(Image[:,1].T, origin = 'lower',
                   cmap = 'inferno',
                   extent=[t[0], t[-1], bins[0], bins[-1]],
                   vmin = Min, vmax = Max,
                   rasterized = True)
ax[1].plot(Time, Means[:,1], c = 'k')
ax[1].set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold')
ax[1].set_ylabel(r'$\mathbf{A.P.L. (\AA^{2})}$', fontsize = 20, fontweight = 'bold')
ax[1].tick_params(labelsize = 18)
forceAspect(ax[1], aspect=1.0)
ax2 = ax[1].twinx()
ax2.plot(t, PR, color = 'w')
ax2.set_ylabel('Pr. (bar)', fontsize = 12, fontweight = 'bold', labelpad = -15)
ax2.set_ylim(pmin, pmax)
ticks = np.linspace(pmin, pmax, 4)
ax2.set_yticks(ticks)
ax2.set_yticklabels(["{:2.0f}".format(int(i)) for i in ticks], fontsize='14')

v1 = np.linspace(Min, Max, 4, endpoint=True)
cb = fig.colorbar(im2, ticks=v1, ax = ax2, shrink = 1, location = 'top')
cb.set_label('P (A.P.L.)',
                fontsize = 11,
                fontweight = 'bold',
                rotation = 0,
                labelpad = 0)
cb.ax.set_xticklabels(["{:4.2f}".format(i) for i in v1], fontsize='14')

fig.tight_layout()

plt.savefig('apl.png', dpi = 300)
plt.savefig('apl.svg', dpi = 300)

