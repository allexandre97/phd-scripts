#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 11:29:32 2022

@author: alexandre
"""

import numpy as np
import matplotlib.pyplot as plt
import cv2 as cv
import os
from natsort import natsorted as ns
from scipy.signal import savgol_filter as sf
import matplotlib.gridspec as gsp

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

linestyles = {'P':'-', 'W':'--'}


files = ns([_ for _ in os.listdir('./DNS/')])

l = len(files)

fig = plt.figure(figsize=(9,12))

gs  = gsp.GridSpec(7,3, figure = fig)

axW = fig.add_subplot(gs[:3,:])
axP = fig.add_subplot(gs[3:6,:])
axT = fig.add_subplot(gs[6,:])

axes = {'P':axP, 'W':axW}

aspect = 0.1

cols = int(len(files)/2)
rows = int(np.ceil(aspect*cols))

time_matrix = np.zeros((rows, cols, 4))

t0 = 0
prev = None
for file in files:
    part = file[4:7]
    if part[0] != prev:
        idx=0
    else:
        idx+=1
    ls   = linestyles[part[0]]
    t1 = t0+5
    label = str(t0)+'ns-'+str(t1)+'ns'
    with open('./DNS/'+file, 'r') as f:
        lines = f.readlines()
        r = np.array([float(line.split()[0]) for line in lines[26:]])
        d = np.array([float(line.split()[1]) for line in lines[26:]])
    f.close()
    d /= np.max(d)
    
    norm_P = (PR[int(PR.shape[0]*(idx/cols)),2,2] - np.min(PR[:,2,2]))/(np.max(PR[:,2,2]) - np.min(PR[:,2,2]))
    rgb = plt.cm.bwr(norm_P)  #(rgb[0]/255, rgb[1]/255, rgb[2]/255)
    
    time_matrix[:,idx] = rgb
    
    axes[part[0]].plot(r, d, color = rgb, label = label)
    t0 = t1
    prev = part[0]
    
pmin = np.min(PR[:,2,2])*(1-0.1)
pmax = np.max(PR[:,2,2])*(1+0.1)

a = (t[-1]-t[0])/abs(80)

axT.imshow(time_matrix, extent = [0, t[-1], -40, 40])
#axT.set_ylim(-40, 40)   
axT.set_xticks(np.linspace(0,500,5))
axT.set_yticks([])
axT.set_xlabel('Time (ns)', fontsize = 14, fontweight = 'bold')
axT.set_ylabel('Pres.', fontsize = 14, fontweight = 'bold')
axT.tick_params(axis = 'x', labelsize = 12)
axT.tick_params(axis = 'y', labelsize = 12)

ax2 = axT.twinx()
ax2.plot(t, PR[:,2,2], color = 'k')
#ax2.set_ylim(-40, 40)
ax2.set_yticks([])

axes['W'].set_title('Solvent', fontsize = 22, fontweight = 'bold')
axes['W'].set_ylabel('Density', fontsize = 20, fontweight = 'bold')
axes['W'].tick_params(labelsize=18)

axes['P'].set_title('PO4 Groups', fontsize = 22, fontweight = 'bold')
axes['P'].set_xlabel('Z (nm)', fontsize = 20, fontweight = 'bold')
axes['P'].set_ylabel('Density', fontsize = 20, fontweight = 'bold')
axes['P'].tick_params(labelsize=18)

fig.tight_layout()