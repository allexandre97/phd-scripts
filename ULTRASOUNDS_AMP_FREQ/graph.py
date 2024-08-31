# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 11:51:04 2023

@author: alexandre
"""

import numpy             as np
import pickle            as pk
import matplotlib.pyplot as plt

from scipy.integrate import simps

with open('./PRES.pk', 'rb') as f:
    PRES = pk.load(f)
f.close()

with open('./FREQ.pk', 'rb') as f:
    FREQ = pk.load(f)
f.close()

with open('./DATA.pk', 'rb') as f:
    DATA = pk.load(f)
f.close()

def GraphData(data    : np.array,
              name    : str,
              outfile : str) -> None:
    
    mean = np.log(data[:,:,0])
    
    stdv = np.log(data[:,:,1] / data[:,:,0])
    
    fig, ax = plt.subplots(1, 2,  figsize = (12, 7),
                           sharex = True, sharey = True)
    
    im1 = ax[0].pcolormesh(PRES, FREQ, mean,
                           cmap = 'inferno')#,
                           #vmin = -25, vmax = 0)
    twinax1 = ax[0].twinx()
    twinax1.set_yticks([])
    
    im2 = ax[1].pcolormesh(PRES, FREQ, stdv,
                           cmap = 'inferno')
    twinax2 = ax[1].twinx()
    twinax2.set_yticks([])
    
    v1 = np.linspace(mean.min(),
                     mean.max(),
                     4, endpoint=True)
    cb1 = fig.colorbar(im1, ticks = v1, ax = twinax1,
                      shrink = 1, location = 'top')
    cb1.set_label(rf'$\mathbf{{ln(\mu_{{{name}}})}}$',
                  fontsize = 18,
                  fontweight = 'bold',
                  rotation = 0,
                  labelpad = 0)
    cb1.ax.set_xticklabels(["{:4.2f}".format(i) for i in v1], fontsize='14')
    
    v2 = np.linspace(stdv.min(),
                     stdv.max(),
                     4, endpoint=True)
    cb2 = fig.colorbar(im2, ticks = v2, ax = twinax2,
                       shrink = 1, location = 'top')
    cb2.set_label(rf'$\mathbf{{ln(\frac{{\sigma_{{{name}}}}}{{\mu_{{{name}}}}})}}$',
                  fontsize = 18,
                  fontweight = 'bold',
                  rotation = 0,
                  labelpad = 0)
    cb2.ax.set_xticklabels(["{:4.2f}".format(i) for i in v2], fontsize='14')
    
    ax[0].set_xlabel('Overpressure / bar',
                     fontsize   = 20,
                      fontweight = 'bold')
    
    ax[0].set_ylabel('Frequency / MHz',
                     fontsize   = 20,
                     fontweight = 'bold')
    
    ax[0].set_xticks(np.arange(5, 55, 5))
    ax[0].set_yticks(np.arange(5, 55, 5))
    
    ax[0].tick_params(labelsize = 18)
    
    
    ax[1].set_xlabel('Overpressure / bar',
                     fontsize   = 20,
                     fontweight = 'bold')
    
    ax[1].set_xticks(np.arange(5, 55, 5))
    
    ax[1].tick_params(labelsize = 18)
    
    fig.tight_layout()
    
    plt.savefig(outfile,
                dpi = 300)
    plt.close()

if __name__ == '__main__':
    
    GraphData(DATA[:,:,0], 'APL',
              'response_apl.svg')
    GraphData(DATA[:,:,1], r'\mathbf{{\kappa}}',
              'response_curv.svg')
    GraphData(DATA[:,:,2], r'Z^{\prime}',
              'response_thick.svg')
    GraphData(DATA[:,:,3], r'S_{1A}',
              'response_ordp_sn1a.svg')
    GraphData(DATA[:,:,4], r'S_{1B}',
              'response_ordp_sn1b.svg')
    GraphData(DATA[:,:,5], r'S_{1C}',
              'response_ordp_sn1c.svg')
    GraphData(DATA[:,:,6], r'S_{2A}',
              'response_ordp_sn2a.svg')
    GraphData(DATA[:,:,7], r'S_{2B}',
              'response_ordp_sn2b.svg')
    GraphData(DATA[:,:,8], r'S_{2C}',
              'response_ordp_sn2c.svg')
    
    names = ['APL', 'CURV', 'THICK', 'SN1A', 'SN1B', 'SN1C', 'SN2A', 'SN2B', 'SN2C']
    
    fig = plt.figure(figsize = (12,9))
    ax = fig.add_subplot(projection = '3d')
    
    for n in range(DATA.shape[2]):
        
        mean = np.log(DATA[:,:,n,0])
        stdv = DATA[:,:,n,1] / DATA[:,:,n,0]
        
        ax.scatter(PRES, FREQ, mean,
                   s = 70*stdv/stdv.max() + 10,
                   label = names[n])
    
    ax.set_xlabel('Pressure / Bar',
                  fontsize   = 20,
                  fontweight = 'bold',
                  labelpad = 10)
    ax.set_ylabel('Frequency / MHz',
                  fontsize   = 20,
                  fontweight = 'bold',
                  labelpad = 10)
    ax.set_zlabel('Response / arb. u.',
                  fontsize   = 20,
                  fontweight = 'bold',
                  labelpad = 20)
    
    ax.tick_params(labelsize = 18)
    
    ax.legend(fontsize = 14)
    
    fig.tight_layout()
    
    
    
    '''def CalcPower(Freq : np.array,
                  Pres : np.array) -> np.array:
        
        POW = np.zeros_like(Freq)
        
        for i in range(Freq.shape[0]):
            for j in range(Freq.shape[1]):
                
                nu = Freq[i,j] * 1E6 * 1E-9
                p  = Pres[i,j]
                
                max_t = 1/nu
                
                time = np.linspace(0, max_t, 10000)
                func = (p * np.sin(2 * np.pi * nu * time))**2
                
                POW[i,j] = simps(func, time) / max_t
                
                
        return POW
                
    PWR = CalcPower(FREQ, PRES)
    
    fig, ax = plt.subplots(figsize = (7,7))
    
    ax.pcolormesh(PRES, FREQ, PWR,
                  cmap = 'inferno')
    
    fig.tight_layout()
    
    RSP = []
    
    for i in range(PWR.shape[0]):
        for j in range(PWR.shape[1]):
            
            tmp = [PWR[i,j]]
            for n in range(DATA.shape[2]):
                tmp.append(DATA[i,j,n,0])
            RSP.append(tmp)
    
    RSP = np.array(RSP)
    
    fig, ax = plt.subplots(figsize = (7, 7))
    
    ax.scatter(RSP[:,0], RSP[:,4])'''
    