#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 13:23:36 2022

@author: alexandre
"""
import sys
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
with open('curv.pickle', 'rb') as f:
    CURV = pickle.load(f)
f.close()

bins = np.linspace(-0.1, 0.1, 100)

'''Timestep parameters to define averaging window (w)
and Time array for mean plotting'''
dt = 20
w  = int(np.ceil(5000/dt))

dt2   = t[-1]/CURV.shape[0]
idxi = int(1000*float(sys.argv[2])/dt2)
idxf = int(1000*float(sys.argv[3])/dt2)

CURV = CURV[idxi:idxf + w]

idxi = int(1000*float(sys.argv[2])/dt1)
idxf = int(1000*float(sys.argv[3])/dt1)

t = t[idxi:idxf]/1000
PR = PR[idxi:idxf]

def Compute(frame, window, CURV):
    
    dummy = CURV[frame:frame+window]
    
    k_up = dummy[:,0].flatten()
    k_dn = dummy[:,1].flatten()
    
    n_up, _ = np.histogram(k_up, bins = bins, density = True)
    n_dn, _ = np.histogram(k_dn, bins = bins, density = True)
    
    return [n_up, n_dn], [[np.mean(k_up), np.std(k_up)],[np.mean(k_dn), np.std(k_dn)]]
    

if __name__ == '__main__':

    import sys    
    import multiprocessing as mp
    from functools import partial
    
    run = partial(Compute,
                  window = w,
                  CURV = CURV)

    Frames = np.arange(CURV.shape[0] - w)

    pool = mp.Pool(int(sys.argv[1]))

    Result = pool.map(run, Frames)

    Khist = np.array([frame[0] for frame in Result])
    K     = np.array([frame[1] for frame in Result])

K     = np.array(K)
Khist = np.array(Khist)

with open('curv_dist.pickle', 'wb') as f:
    pickle.dump(Khist, f)
f.close()

with open('curv_mean.pickle', 'wb') as f:
    pickle.dump(K, f)
f.close()

T = np.linspace(t[0], t[-1], K.shape[0])

Min, Max = np.min([Khist[:,0].min(), Khist[:,1].min()]), np.max([Khist[:,0].max(), Khist[:,1].max()])*(1-0.25)

fig, ax = plt.subplots(1, 2, figsize = (12,5))

ax[0].set_title('Upper Leaflet', fontsize = 22, fontweight = 'bold', pad = 60)
im1 = ax[0].imshow(Khist[:,0].T, 
                   cmap = 'inferno',
                   origin = 'lower',
                   extent = [t[0], t[-1], bins[0], bins[-1]],
                   vmin = Min, vmax = Max)
ax[0].plot(T, K[:,0,0], c = 'k')
ax[0].set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold')
ax[0].set_ylabel('K ($\mathbf{nm^{-1}}$)', fontsize = 20, fontweight = 'bold')
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
cb = fig.colorbar(im1, ticks=v1, ax = ax2, shrink = 1, location ='top')
cb.set_label(r'$\mathbf{P(K)}$',
                fontsize = 11,
                fontweight = 'bold',
                rotation = 0,
                labelpad = -15)
cb.ax.set_xticklabels(["{:4.2f}".format(i) for i in v1], fontsize='14')

ax[0].set_xlim(float(sys.argv[2]), float(sys.argv[3]))


ax[1].set_title('Lower Leaflet', fontsize = 22, fontweight = 'bold', pad = 60)
im2 = ax[1].imshow(Khist[:,1].T, 
                   cmap = 'inferno',
                   origin = 'lower',
                   extent = [t[0], t[-1], bins[0], bins[-1]],
                   vmin = Min, vmax = Max)
ax[1].plot(T, K[:,1,0], c = 'k')
ax[1].set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold')
ax[1].set_ylabel('K ($\mathbf{nm^{-1}}$)', fontsize = 20, fontweight = 'bold')
ax[1].tick_params(labelsize = 18)
#ax[1].set_aspect(a)
ax2 = ax[1].twinx()
ax2.plot(t, PR, color = 'w')
ax2.set_ylabel('Pr. (bar)', fontsize = 12, fontweight = 'bold', labelpad = -15)
ax2.set_ylim(pmin, pmax)
ticks = np.linspace(pmin, pmax, 4)
ax2.set_yticks(ticks)
ax2.set_yticklabels(["{:2.0f}".format(int(i)) for i in ticks], fontsize='14')
forceAspect(ax[1], aspect=1.0)

v2 = np.linspace(Min, Max, 4, endpoint=True)
cb2 = fig.colorbar(im2, ticks=v2, ax = ax2, shrink = 1, location = 'top')
cb2.set_label(r'$\mathbf{P(K)}$',
                fontsize = 11,
                fontweight = 'bold',
                rotation = 0,
                labelpad = -15)
cb2.ax.set_xticklabels(["{:4.2f}".format(i) for i in v1], fontsize='14')

ax[1].set_xlim(float(sys.argv[2]), float(sys.argv[3]))

fig.tight_layout()

plt.savefig('curv.png', dpi = 300)
plt.savefig('curv.svg', dpi = 300)

