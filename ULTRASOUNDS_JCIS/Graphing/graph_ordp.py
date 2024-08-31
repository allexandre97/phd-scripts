#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 13:16:52 2022

@author: alexandre
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsp
import pickle
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
with open(f'ordpars_{sys.argv[1]}.pickle', 'rb') as f:
    ORD = pickle.load(f)
f.close()

'''Timestep parameters to define averaging window (w)
and Time array for mean plotting'''
dt = 20
w  = int(np.ceil(5000/dt))

dt2   = t[-1]/ORD.shape[0]
idxi = int(1000*float(sys.argv[3])/dt2)
idxf = int(1000*float(sys.argv[4])/dt2)

ORD = ORD[idxi:idxf + w]

idxi = int(1000*float(sys.argv[3])/dt1)
idxf = int(1000*float(sys.argv[4])/dt1)

t = t[idxi:idxf]/1000
PR = PR[idxi:idxf]

def Compute(frame, window, ORD):
    
    data1_up = ORD[frame:frame+window,0,0]
    data1_dn = ORD[frame:frame+window,0,1]

    data2_up = ORD[frame:frame+window,1,0]
    data2_dn = ORD[frame:frame+window,1,1]
    
    mean1_up = np.mean(data1_up, axis = 0)
    mean1_dn = np.mean(data1_dn, axis = 0)
    
    mean1_up = np.mean(mean1_up, axis = 0)
    mean1_dn = np.mean(mean1_dn, axis = 0)
    
    ord1_up = 0.5*(3*mean1_up - 1)
    ord1_dn = 0.5*(3*mean1_dn - 1)

    mean2_up = np.mean(data2_up, axis = 0)
    mean2_dn = np.mean(data2_dn, axis = 0)
    
    mean2_up = np.mean(mean2_up, axis = 0)
    mean2_dn = np.mean(mean2_dn, axis = 0)
    
    ord2_up = 0.5*(3*mean2_up - 1)
    ord2_dn = 0.5*(3*mean2_dn - 1)
    
    return [[ord1_up, ord1_dn],[ord2_up, ord2_dn]]

if __name__ == '__main__':
    
    import multiprocessing as mp
    from functools import partial
    
    run = partial(Compute,
                  window = w,
                  ORD = ORD)

    pool = mp.Pool(int(sys.argv[2]))

    Frames = np.arange(ORD.shape[0] - w)

    Result = pool.map(run, Frames)

    OrdPars = np.array(Result)


Time = np.linspace(t[0], t[-1], OrdPars.shape[0])

fig = plt.figure(figsize = (12,9))

gs  = gsp.GridSpec(9,8, figure = fig)

ax1UP = fig.add_subplot(gs[0:4,0:4])
ax1DN = fig.add_subplot(gs[0:4,4:])

ax2UP = fig.add_subplot(gs[4:-1,0:4])
ax2DN = fig.add_subplot(gs[4:-1,4:])

axPUP  = fig.add_subplot(gs[-1,0:4])
axPDN  = fig.add_subplot(gs[-1,4:])

ax1UP.set_title('Upper Leaflet', fontweight = 'bold', fontsize = 22)
ax1DN.set_title('Lower Leaflet', fontweight = 'bold', fontsize = 22)

aux1 = ax1DN.twinx()
aux1.set_ylabel('SN1\nTail', fontsize = 20, fontweight = 'bold', rotation = 0, labelpad = 35)
aux1.set_yticks([])

aux2 = ax2DN.twinx()
aux2.set_ylabel('SN2\nTail', fontsize = 20, fontweight = 'bold', rotation = 0, labelpad = 35)
aux2.set_yticks([])

ax1UP.plot(Time, OrdPars[:,0,0,0], c = 'royalblue')
ax1UP.plot(Time, OrdPars[:,0,0,1], c = 'firebrick')
ax1UP.plot(Time, OrdPars[:,0,0,2], c = 'seagreen')

ax1UP.set_ylabel(r'$\mathbf{S_{seg}}$', fontsize = 20, fontweight = 'bold')
ax1UP.set_xticklabels([])
ax1UP.set_xlim(t[0], t[-1])
ax1UP.tick_params(labelsize = 18)

ax2UP.plot(Time, OrdPars[:,1,0,0], c = 'royalblue')
ax2UP.plot(Time, OrdPars[:,1,0,1], c = 'firebrick')
ax2UP.plot(Time, OrdPars[:,1,0,2], c = 'seagreen')

ax2UP.set_ylabel(r'$\mathbf{S_{seg}}$', fontsize = 20, fontweight = 'bold')
ax2UP.set_xticklabels([])
ax2UP.set_xlim(t[0], t[-1])
ax2UP.tick_params(labelsize = 18)

axPUP.plot(t, PR, color = 'k')
axPUP.set_ylabel('Pres.', fontsize = 12, fontweight = 'bold', labelpad = 0)
axPUP.set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold', labelpad = 0)
axPUP.set_ylim(pmin, pmax)
ticks = np.linspace(pmin, pmax, 2)
axPUP.set_yticks([])
axPUP.set_yticklabels([])
axPUP.tick_params(labelsize = 20)
axPUP.set_xlim(t[0], t[-1])

ax1DN.plot(Time, OrdPars[:,0,1,0], c = 'royalblue')
ax1DN.plot(Time, OrdPars[:,0,1,1], c = 'firebrick')
ax1DN.plot(Time, OrdPars[:,0,1,2], c = 'seagreen')

ax1DN.set_xlim(t[0], t[-1])
ax1DN.set_xticklabels([])
ax1DN.set_yticklabels([])

ax2DN.plot(Time, OrdPars[:,1,1,0], c = 'royalblue')
ax2DN.plot(Time, OrdPars[:,1,1,1], c = 'firebrick')
ax2DN.plot(Time, OrdPars[:,1,1,2], c = 'seagreen')

ax2DN.set_xlim(t[0], t[-1])
ax2DN.set_xticklabels([])
ax2DN.set_yticklabels([])

axPDN.plot(t, PR, color = 'k')
axPDN.set_xlabel('Time (ns)', fontsize = 20, fontweight = 'bold', labelpad = 0)
axPDN.set_ylim(pmin, pmax)
ticks = np.linspace(pmin, pmax, 2)
axPDN.set_yticks([])
axPDN.set_yticklabels([])
axPDN.tick_params(labelsize = 20)
axPDN.set_xlim(t[0], t[-1])

ax1UP.set_ylim(0, 1)
ax2UP.set_ylim(0, 1)
ax1DN.set_ylim(0, 1)
ax2DN.set_ylim(0, 1)


fig.tight_layout()

plt.savefig(f'ordp_{sys.argv[1]}.png', dpi = 300)
plt.savefig(f'ordp_{sys.argv[1]}.svg', dpi = 300)
