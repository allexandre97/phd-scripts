#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 13:23:36 2022

@author: alexandre
"""

import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter as sf

def forceAspect(ax,aspect):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
    
'''Load pressure tensor from properties xvg file'''
with open('pressure.xvg', 'r') as f:
    TP = np.loadtxt(f, comments = ('#', '@')
f.close()

'''Store pressure tensor in one array'''
PRE = TP[:,1]
t   = TP[:,0]

PR = sf(PRE, 80*50+1, 2, axis = 0)

pmin = np.min(PR[:,2,2])*(1-0.1)
pmax = np.max(PR[:,2,2])*(1+0.1)

with open('thickness.pickle', 'rb') as f:
    THI = pickle.load(f)*0.1
f.close()

with open('gridX.pickle', 'rb') as f:
    gX = pickle.load(f)
f.close()

with open('gridY.pickle', 'rb') as f:
    gY = pickle.load(f)
f.close()

dt = 20
w  = int(np.ceil(5000/dt))

bns = np.linspace(2,7,100)

def Compute(Frame, window, THI):
    
    dummy = THI[Frame:Frame+window]
    
    thick = dummy.flatten()
    
    n, _ = np.histogram(thick, bins = bns, density = True)
    
    return n, [np.mean(thick), np.std(thick)]


if __name__ == '__main__':
    
    import sys
    import os
    import multiprocessing as mp
    from functools import partial
    
    chunksize = sys.argv[1]*50

    maxiter   = int(np.ceil((THI.shape[0]-w)/chunksize))

    Result = []

    with mp.Pool(sys.argv[1]) as pool:

        for iteration in range(maxiter):

            begin = iteration*chunksize

            if (iteration+1)*chunksize < THI.shape[0]-w:
                end = (iteration+1)*chunksize
            else:
                end = -1

            Frames = np.arange(THI[begin:end].shape[0])

            run = partial(Compute,
                          window = w,
                          THI = THI[begin:end])

            result = pool.map(run, Frames)

            for frame in result:
                Result.append(frame)

    pool.close()
    pool.join()
    Thist     = np.array([frame[0] for frame in Result])
    T = np.array([frame[1] for frame in Result])    
    

with open('thickness_dist.pickle', 'wb') as f:
    pickle.dump(Thist, f)
f.close()

with open('thickness_mean.pickle', 'wb') as f:
    pickle.dump(T, f)
f.close()

Time = 0.02*(np.linspace(0, T.shape[0], T.shape[0]))

a = (Time[-1]-Time[0])/(bns[-1]-bns[0])

fig, ax = plt.subplots(figsize = (9,9))

im1 = ax.imshow(Thist.T, cmap = 'inferno', origin = 'lower', extent = [t[0], t[-1], bns[0], bns[-1]])
ax.plot(Time, T[:,0], color = 'k')
#ax[0].plot(Time, T[:,0,0], c = 'firebrick')
#ax[0].fill_between(Time, T[:,0,0]-T[:,0,1], T[:,0,0]+T[:,0,1], color = 'firebrick', alpha = 0.2)
ax.set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold')
ax.set_ylabel('Th. ($\mathbf{nm}$)', fontsize = 20, fontweight = 'bold')
ax.tick_params(labelsize = 18)
forceAspect(ax, aspect=1.0)

ax2 = ax.twinx()
ax2.plot(t, PR[:,2,2], color = 'w')
ax2.set_ylabel('Pr. (bar)', fontsize = 12, fontweight = 'bold', labelpad = -15)
ax2.set_ylim(pmin-10, pmax)
ticks = np.linspace(pmin-10, pmax, 4)
ax2.set_yticks(ticks)
ax2.set_yticklabels(["{:2.0f}".format(int(i)) for i in ticks], fontsize='14')

v1 = np.linspace(Thist.min(), Thist.max(), 4, endpoint=True)
cb = fig.colorbar(im1, ticks=v1, ax = ax2, shrink = 1, location = 'top')
cb.set_label(r'$\mathbf{P(Th.)}$',
                fontsize = 11,
                fontweight = 'bold',
                rotation = 0,
                labelpad = -15)
cb.ax.set_xticklabels(["{:4.2f}".format(i) for i in v1], fontsize='14')

ax.set_xlim(sys.argv[2], sys.argv[3])

fig.tight_layout()

