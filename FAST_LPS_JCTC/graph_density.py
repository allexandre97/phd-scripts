# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 12:11:32 2023

@author: alexandre
"""

import numpy             as np
import pickle            as pk
import matplotlib.pyplot as plt

class Dens:

    def __init__(self, BINS, DENSITY):

        self.BINS    = BINS
        self.DENSITY = DENSITY

FULLDATA = []

for replica in range(1,4,1):
    
    REP = []
    
    for folder in ['NoFast', 'HalfFast', 'AllFast']:
        
        fullpath = f'Rep{replica}/{folder}'

        with open(f'{fullpath}/density.pk', 'rb') as f:
            DENS = pk.load(f)
        f.close()


        dt = 5000000/DENS.DENSITY.shape[0]
        W  = int(np.ceil(50000/dt))
        
        idx_0 = int(np.ceil(0/dt))
        idx_1 = int(np.ceil(500000/dt)) - W
        
        skip = 10
        
        WINDOWS = np.lib.stride_tricks.sliding_window_view(DENS.DENSITY, W, axis = 0)
        
        colors = ['#648fff', '#fe6100', '#000000', '#dc267f']
        
        fig, ax = plt.subplots(figsize = (7,7))
        
        for t in range(idx_0, idx_1, skip):
        
            if t == idx_1-skip:
        
                MEAN_DENS = np.mean(DENS.DENSITY[idx_0:idx_1], axis = 0)
                
                print(MEAN_DENS.shape)
                
                REP.append([DENS.BINS, MEAN_DENS])
        
                X = DENS.BINS
        
                alpha = 1
        
                ax.plot(X, MEAN_DENS[0],
                        c = colors[0], ls = '-', label = 'Polymyxin', alpha = alpha)
        
                if (MEAN_DENS[1] != 0).any():
                    ax.plot(X, MEAN_DENS[1],
                            c = colors[1], ls = '-', label = 'REMP-Head', alpha = alpha)
        
                    ax.plot(X, MEAN_DENS[2],
                            c = colors[2], ls = '-', label = 'REMP-Tails', alpha = alpha)
        
                if (MEAN_DENS[3] != 0).any():
                    ax.plot(X, MEAN_DENS[3],
                            c = colors[1], ls = '--', label = 'Fast REMP-Head', alpha = alpha)
        
                    ax.plot(X, MEAN_DENS[4],
                            c = colors[2], ls = '--', label = 'Fast REMP-Tails', alpha = alpha)
        
                ax.plot(X, MEAN_DENS[5],
                        c = colors[3], ls = '-', label = 'Lipids', alpha = alpha)
        
            else:
        
        
                MEAN_DENS = np.mean(WINDOWS[t], axis = -1)
        
                X = DENS.BINS
        
                alpha = 0.01
        
                ax.plot(X, MEAN_DENS[0],
                        c = colors[0], alpha = alpha)
        
                ax.plot(X, MEAN_DENS[1],
                        c = colors[1], alpha = alpha)
        
                ax.plot(X, MEAN_DENS[2],
                        c = colors[2], alpha = alpha)
        
                ax.plot(X, MEAN_DENS[3],
                        c = colors[1], ls = '--', alpha = alpha)
        
                ax.plot(X, MEAN_DENS[4],
                        c = colors[2], ls = '--', alpha = alpha)
        
                ax.plot(X, MEAN_DENS[5],
                        c = colors[3], alpha = alpha)
        
        
        ax.set_ylim(0, 3.5)
        ax.set_xlabel(r'Z / $\mathbf{\AA}$',
                      fontsize = 20, fontweight = 'bold')
        ax.set_ylabel(r'$\mathbf{\rho}$ / $\mathbf{\frac{amu}{\AA^3}}$',
                      fontsize = 20, fontweight = 'bold')
        ax.tick_params(labelsize = 18)
        
        ax.legend(fontsize = 14, loc = 'upper right')
        
        fig.tight_layout()
        
        plt.savefig(f'{fullpath}/density_initial.svg',
                    format = 'svg',
                    dpi = 600)
    
    FULLDATA.append(REP)

from scipy.interpolate import interp1d

def InterpolateFulldata(FULLDATA, fast):
    
    names = {0:'NoFast', 1:'HalfFast', 2:'AllFast'}
    
    err = 1
    
    X = np.unique(np.concatenate((FULLDATA[0][fast][0], FULLDATA[1][fast][0], FULLDATA[2][fast][0])))
    
    Y1 = interp1d(FULLDATA[0][fast][0], FULLDATA[0][fast][1],
                  kind = 'linear', fill_value = 'extrapolate')(X)
    Y2 = interp1d(FULLDATA[1][fast][0], FULLDATA[1][fast][1],
                  kind = 'linear', fill_value = 'extrapolate')(X)
    Y3 = interp1d(FULLDATA[2][fast][0], FULLDATA[2][fast][1],
                  kind = 'linear', fill_value = 'extrapolate')(X)
    
    Y = np.array([np.mean([Y1, Y2, Y3], axis = 0),
                  np.std([Y1, Y2, Y3], axis = 0)])
    
    colors = ['#648fff', '#fe6100', '#000000', '#dc267f']
    
    fig, ax = plt.subplots(figsize = (7,7))

    alpha = 1

    ax.plot(X, Y[0,0],
            c = colors[0], ls = '-', label = 'Polymyxin', alpha = alpha)
    ax.fill_between(X, Y[0,0] - err*Y[1,0], Y[0,0] + err*Y[1,0],
                    color = colors[0], edgecolor = None, alpha = 0.25)

    if (Y[0,1] != 0).any():
        ax.plot(X, Y[0,1],
                c = colors[1], ls = '-', label = 'REMP-Head', alpha = alpha)
        ax.fill_between(X, Y[0,1] - err*Y[1,1], Y[0,1] + err*Y[1,1],
                        color = colors[1],  edgecolor = None, alpha = 0.25)
        
        ax.plot(X, Y[0,2],
                c = colors[2], ls = '-', label = 'REMP-Tails', alpha = alpha)
        ax.fill_between(X, Y[0,2] - err*Y[1,2], Y[0,2] + err*Y[1,2],
                        color = colors[2],  edgecolor = None, alpha = 0.25)
        
    if (Y[0,3] != 0).any():
        ax.plot(X, Y[0,3],
                c = colors[1], ls = '--', label = 'Fast REMP-Head', alpha = alpha)
        ax.fill_between(X, Y[0,3] - err*Y[1,3], Y[0,3] + err*Y[1,3],
                        color = colors[1],  edgecolor = None, alpha = 0.25)
        
        ax.plot(X, Y[0,4],
                c = colors[2], ls = '--', label = 'Fast REMP-Tails', alpha = alpha)
        ax.fill_between(X, Y[0,4] - err*Y[1,4], Y[0,4] + err*Y[1,4],
                        color = colors[2],  edgecolor = None, alpha = 0.25)
        
    ax.plot(X, Y[0,5],
            c = colors[3], ls = '-', label = 'Lipids', alpha = alpha)
    ax.fill_between(X, Y[0,5] - err*Y[1,5], Y[0,5] + err*Y[1,5],
                    color = colors[3],  edgecolor = None, alpha = 0.25)
    
    ax.set_ylim(0, 3.5)
    ax.set_xlabel(r'Z / $\mathbf{\AA}$',
                  fontsize = 20, fontweight = 'bold')
    ax.set_ylabel(r'$\mathbf{\rho}$ / $\mathbf{\frac{amu}{\AA^3}}$',
                  fontsize = 20, fontweight = 'bold')
    ax.tick_params(labelsize = 18)
    
    ax.legend(fontsize = 14, loc = 'upper right')
    
    fig.tight_layout()
    
    plt.savefig(f'avg_density_initial_{names[fast]}.svg',
                format = 'svg',
                dpi = 600)

InterpolateFulldata(FULLDATA, 0)
InterpolateFulldata(FULLDATA, 1)
InterpolateFulldata(FULLDATA, 2)


