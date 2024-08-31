#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 10:45:56 2022

@author: alexandre
"""

import sys
import pickle 
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter as sf

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
with open('flipflop.pickle', 'rb') as f:
    FL = pickle.load(f)
f.close()

dt2   = t[-1]/FL.shape[0]
idxi = int(1000*float(sys.argv[1])/dt2)
idxf = int(1000*float(sys.argv[2])/dt2)

FL = FL[idxi:idxf]

idxi = int(1000*float(sys.argv[1])/dt1)
idxf = int(1000*float(sys.argv[2])/dt1)

t = t[idxi:idxf]/1000
PR = PR[idxi:idxf]

T = np.linspace(t[0], t[-1],FL.shape[0])

import matplotlib.cm as cm

fig, ax = plt.subplots(1,2,figsize = (12,5), sharex = True)

L = FL.shape[2]

for l in range(L):
    
    ax[0].scatter(T, 0.1*FL[:,0,l], color = cm.inferno(l/L), alpha = 0.4, rasterized = True)

ax2 = ax[0].twinx()
ax2.plot(t, PR, color = 'firebrick')
ax2.set_ylabel('Pr. (bar)', fontsize = 12, fontweight = 'bold', labelpad = -15)
ax2.set_ylim(pmin, pmax)
ticks = np.linspace(pmin, pmax, 4)
ax2.set_yticks(ticks)
ax2.set_yticklabels(["{:2.0f}".format(int(i)) for i in ticks], fontsize='14')


for l in range(L):
    
    ax[1].scatter(T, 0.1*FL[:,1,l], color = cm.inferno(l/L), alpha = 0.4, rasterized = True)

ax2 = ax[1].twinx()
ax2.plot(t, PR, color = 'firebrick')
ax2.set_ylabel('Pr. (bar)', fontsize = 12, fontweight = 'bold', labelpad = -15)
ax2.set_ylim(pmin, pmax)
ticks = np.linspace(pmin, pmax, 4)
ax2.set_yticks(ticks)
ax2.set_yticklabels(["{:2.0f}".format(int(i)) for i in ticks], fontsize='14')

ax[0].set_title('Upper Leaflet', fontsize = 22, fontweight = 'bold')
ax[1].set_title('Lower Leaflet', fontsize = 22, fontweight = 'bold')
ax[0].set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold')
ax[0].set_ylabel('Distance to C.O.M. (nm)', fontsize = 20, fontweight = 'bold')
ax[1].set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold')
ax[1].set_ylabel('Distance to C.O.M. (nm)', fontsize = 20, fontweight = 'bold')
ax[0].set_xticks(np.linspace(0,500,4))
ax[1].set_xticks(np.linspace(0,500,4))
ax[0].tick_params(labelsize = 20)
ax[1].tick_params(labelsize = 20)

ax[0].set_xlim(float(sys.argv[1]), float(sys.argv[2]))
ax[1].set_xlim(float(sys.argv[1]), float(sys.argv[2]))


fig.tight_layout()
plt.savefig('flipflop.svg', dpi = 300)
plt.savefig('flipflop.png', dpi = 300)
