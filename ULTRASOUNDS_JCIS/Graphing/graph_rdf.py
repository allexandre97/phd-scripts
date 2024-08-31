#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 16:21:37 2022

@author: alexandre
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsp
import cv2 as cv
from natsort import natsorted
from scipy.signal import savgol_filter as sf


def forceAspect(ax,aspect):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

    
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

pmin = np.min(PR[:,2,2])*(1+0.1)
pmax = np.max(PR[:,2,2])*(1+0.1)

# Get all the filenames
files = natsorted([_ for _ in os.listdir('./RDF/') if _.startswith('rdf') and _.endswith('.xvg')])

# Get all the studied lipids
parts = []
for fil in files:
    p1 = fil[4:6]
    p2 = fil[7:9]
    if not p1 in parts:
        parts.append(p1)
    if not p2 in parts:
        parts.append(p2)
parts = natsorted(parts,reverse = True)

dt = 1.25

fig  = plt.figure(figsize = (9,9))
fig2 = plt.figure(figsize = (9,9))


l = len(parts)
d = 4
s = d*l

gs  = gsp.GridSpec(s+int(s/4),s, figure = fig)
gs2 = gsp.GridSpec(s,s, figure = fig2)

aspect = 0.1

cols = int(len(files)/len(parts))
rows = int(np.ceil(aspect*cols))

time_matrix = np.zeros((rows, cols, 4))


row, col = 0, 0
for i in range(len(parts)):
    for j in range(i, len(parts)):
        
        ps = natsorted([parts[i], parts[j]])
        par = ps[0]+'-'+ps[1]
        ax  = fig.add_subplot(gs[row*d:d*(row+1),col*d:d*(col+1)])
        ax2 = fig2.add_subplot(gs2[row*d:d*(row+1),col*d:d*(col+1)])
        
        kk = []
        
        idx = 0
        
        for file in files:
            if file[4:9] == par:
                                
                with open('./RDF/'+file, 'r') as f:
                    lines = f.readlines()
                    d   = [float(line.split()[0]) for n, line in enumerate(lines) if n > 26 and (0.35 <= float(line.split()[0]) <= 1.5) ]
                    rdf = [float(line.split()[1]) for n, line in enumerate(lines) if n > 26 and (0.35 <= float(line.split()[0]) <= 1.5) ]
                    d   = np.array(d)
                    rdf = np.array(rdf)
                f.close()
                
                rdf = sf(rdf, 51, 3)
                kk.append(rdf)
                '''
                hsv = np.uint8([[[90 + 120*(idx/cols), 230, 200]]])
                rgb = cv.cvtColor(hsv, cv.COLOR_HSV2RGB)[0,0]
                '''
                
                norm_P = (PR[int(PR.shape[0]*(idx/cols)),2,2] - np.min(PR[:,2,2]))/(np.max(PR[:,2,2]) - np.min(PR[:,2,2]))
                
                rgb = plt.cm.bwr(norm_P)  #(rgb[0]/255, rgb[1]/255, rgb[2]/255)

                
                time_matrix[:,idx] = rgb
                
                label = str(dt*(idx))+'ns - '+str(dt+dt*(idx))+'ns'
                
                ax.plot(d, rdf, color = rgb, label = label)
                
                idx += 1
        
        ax.set_title(par, fontsize = 22/l, fontweight = 'bold')
        ax.set_xlabel('R (nm)', fontsize = 20/l, fontweight = 'bold')
        ax.set_ylabel('P(R)', fontsize = 20/l, fontweight = 'bold')
        ax.tick_params(labelsize=18/l)
        #ax.legend()
        
        kk = np.array(kk)
        
        a = 100/d[-1]
        
        ax2.set_title(par, fontsize = 22/l, fontweight = 'bold', pad = 90)
        im1 = ax2.imshow(kk.T, origin = 'lower',
                   cmap = 'inferno',
                   extent = [0, 500, d[0], d[-1]])
        ax2.set_xlabel('Time (ns)', fontsize = 20/l, fontweight = 'bold')
        ax2.set_ylabel('R (nm)', fontsize = 20/l, fontweight = 'bold')
        ax2.tick_params(labelsize = 18/l)
        axP = ax2.twinx()
        axP.plot(t, PR[:,2,2], color = 'w')
        axP.set_ylabel('P (bar)', fontsize = 12, fontweight = 'bold', labelpad = -15)
        axP.set_ylim(pmin, pmax)
        ticks = np.linspace(pmin, pmax, 4)
        axP.set_yticks(ticks)
        axP.set_yticklabels(["{:2.0f}".format(int(i)) for i in ticks], fontsize='14')
        forceAspect(ax2, aspect=1.0)

        v1 = np.linspace(kk.min(), kk.max(), 4, endpoint=True)
        cb = fig.colorbar(im1, ticks=v1, ax = axP, shrink = 1, location = 'top')
        cb.set_label(r'$\mathbf{P(R)}$',
                        fontsize = 11,
                        fontweight = 'bold',
                        rotation = 0,
                        labelpad = -15)
        cb.ax.set_xticklabels(["{:4.2f}".format(i) for i in v1], fontsize='14')
        
        if col < (l-1):
            col += 1
        else:
            row += 1
            col = row

axT = fig.add_subplot(gs[-1:,:])
axT.imshow(time_matrix, extent = [0, t[-1], -40, 40])
#axT.set_ylim(-40, 40)   
axT.set_xticks(np.linspace(0,500,5))
axT.set_yticks([])
axT.set_xlabel('Time (ns)', fontsize = 14, fontweight = 'bold')
axT.set_ylabel('Pres.', fontsize = 14, fontweight = 'bold')
axT.tick_params(axis = 'x', labelsize = 12)
axT.tick_params(axis = 'y', labelsize = 12)
axT.set_aspect(2)
axP = axT.twinx()
axP.plot(t, PR[:,2,2], color = 'k')
#ax2.set_ylim(-40, 40)
axP.set_yticks([])

ax2.set_ylim(d[0],d[-1])

fig.tight_layout()
fig2.tight_layout()

fig.savefig('rdf_1.png', dpi = 300)
fig.savefig('rdf_1.svg', dpi = 300)

fig2.savefig('rdf_2.png', dpi = 300)
fig2.savefig('rdf_2.svg', dpi = 300)
