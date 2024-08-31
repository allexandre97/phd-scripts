#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 10:11:40 2024

@author: alexandre
"""

import numpy             as np
import pickle            as pk
import matplotlib.pyplot as plt

from scipy.stats         import t, ttest_ind_from_stats
from matplotlib.gridspec import GridSpec as GS

sig_level = 0.001

MEMBRANES = ['POPC',
             'POPE',
             'POPG',
             'POPS']


def Graph(P : str = 'APL') -> None:
    
    PROP = {'APL':0,
            'CUR':1,
            'THI':2,
            'SN1A':3,
            'SN1B':4,
            'SN1C':5,
            'SN2A':6,
            'SN2B':7,
            'SN2C':8}
    
    PROPNAME = {'APL':'APL',
                'CUR':r'\kappa',
                'THI':r'Z^{\prime}',
                'SN1A':r'S_{1A}',
                'SN1B':r'S_{1B}',
                'SN1C':r'S_{1C}',
                'SN2A':r'S_{2A}',
                'SN2B':r'S_{2B}',
                'SN2C':r'S_{2C}'}
    
    fig = plt.figure(figsize = (10, 10))
    
    gs = GS(5, 4, figure = fig, height_ratios = (1,1,1,1,0.25))
    
    ticks = np.linspace(10, 50, 5)
    
    n = 0
    for i in range(4):
        
        for j in range(4):
            
            data = DIFF[n]
            pval = PVAL[n]
            
            mask = np.zeros((10, 10, 4))
            
            mask[pval[:,:,PROP[P.upper()],0] > sig_level, 3] = 0.5
            
            # if i==0 and j==0:
                
            ax = fig.add_subplot(gs[i:i+1,j:j+1])
            
            # else:
                
            #     ax = fig.add_subplot(gs[i:i+1,j:j+1], sharex = axes[0], sharey = axes[0])
            
            if j == 0:
                
                ax.set_ylabel(r'$\mathbf{\nu}$ / MHz',
                              fontsize = 12, fontweight = 'bold')
                ax.set_yticks(ticks)
            
            elif j == 3:
                
                twin = ax.twinx()
                twin.set_ylabel(f'{MEMBRANES[i]}',
                                fontsize = 16, fontweight = 'bold',
                                rotation = 0, labelpad = 27)
                twin.set_yticks([])
                ax.set_yticks([])
                
            else:
                ax.set_yticks([])
                
            if i == 3:
                ax.set_xlabel('P / bar',
                              fontsize = 12, fontweight = 'bold')
                ax.set_xticks(ticks)
            
            elif i == 0:
                ax.set_title(f'{MEMBRANES[j]}',
                             fontsize = 16, fontweight = 'bold')
                ax.set_xticks([])
                
            else:
                ax.set_xticks([])
            
                
            minval = MIN[PROP[P.upper()]]
            maxval = MAX[PROP[P.upper()]]
        
                
            img = ax.pcolormesh(PRES, FREQ, data[:,:,PROP[P.upper()],0],
                                cmap = 'bwr',
                                vmin = minval, vmax = maxval)
            mask_img = ax.pcolormesh(PRES, FREQ, mask)
            
            n += 1
    
    
    cbar_ticks = np.linspace(minval, maxval, 4) 
    ax_cbar    = fig.add_subplot(gs[-1,:])
    
    cbar = fig.colorbar(img, cax = ax_cbar,
                        ticks = cbar_ticks,
                        location = 'top')
    cbar.set_label(r'$\mathbf{\Delta \mu_{%s}}$' % (PROPNAME[P.upper()]),
                   fontsize = 16, fontweight = 'bold',
                   rotation = 0, labelpad = -50)
    
    cbar.ax.tick_params(top = False, labeltop = False,
                        bottom = True, labelbottom = True)
    cbar.ax.set_xticklabels(['{:4.2E}'.format(val) for val in cbar_ticks],
                            fontsize = 12)
    
    fig.tight_layout()
    
    plt.savefig(f'diff_{P.upper()}.svg',
                dpi = 600)
#%%
def GraphZoom(P : str = 'APL',
              I : int = 2,
              J : int = 1) -> None:
    
    PROP = {'APL':0,
            'CUR':1,
            'THI':2,
            'SN1A':3,
            'SN1B':4,
            'SN1C':5,
            'SN2A':6,
            'SN2B':7,
            'SN2C':8}
    
    PROPNAME = {'APL':'APL',
                'CUR':r'\kappa',
                'THI':r'Z^{\prime}',
                'SN1A':r'S_{1A}',
                'SN1B':r'S_{1B}',
                'SN1C':r'S_{1C}',
                'SN2A':r'S_{2A}',
                'SN2B':r'S_{2B}',
                'SN2C':r'S_{2C}'}
    
    n = I*4 + J
            
    ticks = np.linspace(10, 50, 5)
    
    fig, ax = plt.subplots(figsize = (7, 7))
    
    data = DIFF[n]
    pval = PVAL[n]
    
    minval = MIN[PROP[P.upper()]]
    maxval = MAX[PROP[P.upper()]]
    
    mask = np.zeros((10, 10, 4))
    
    mask[pval[:,:,PROP[P.upper()],0] > sig_level, 3] = 0.5
    
    img = ax.pcolormesh(PRES, FREQ, data[:,:,PROP[P.upper()],0],
                        cmap = 'bwr',
                        vmin = minval, vmax = maxval)
    mask_img = ax.pcolormesh(PRES, FREQ, mask)
    
    ax.set_ylabel(r'$\mathbf{\nu}$ / MHz',
                  fontsize = 14, fontweight = 'bold')
    ax.set_yticks(ticks)
    
    ax.set_xlabel('P / bar',
                  fontsize = 14, fontweight = 'bold')
    ax.set_xticks(ticks)
    
    ax.tick_params(labelsize = 12)
    
    
    ax.set_title(f'{MEMBRANES[J]}',
                 fontsize = 16, fontweight = 'bold',
                 pad = 85)
    twin = ax.twinx()
    twin.set_ylabel(f'{MEMBRANES[I]}',
                    fontsize = 16, fontweight = 'bold',
                    rotation = 0, labelpad = 27)
    twin.set_yticks([])
    
    
    
    cbar_ticks = np.linspace(minval, maxval, 4) 
    
    cbar = fig.colorbar(img, ax = twin,
                        ticks = cbar_ticks,
                        location = 'top')
    cbar.set_label(r'$\mathbf{\Delta \mu_{%s}}$' % (PROPNAME[P.upper()]),
                   fontsize = 16, fontweight = 'bold',
                   rotation = 0, labelpad = 0)
    
    cbar.ax.tick_params(top    = True,  labeltop    = True,
                        bottom = False, labelbottom = False)
    cbar.ax.set_xticklabels(['{:4.2E}'.format(val) for val in cbar_ticks],
                            fontsize = 12)
    
    fig.tight_layout()
    
    # plt.savefig(f'diff_{P.upper()}.svg',
    #             dpi = 600)
#%%

def GraphPV(P : str = 'APL') -> None:
    
    PROP = {'APL':0,
            'CUR':1,
            'THI':2,
            'SN1A':3,
            'SN1B':4,
            'SN1C':5,
            'SN2A':6,
            'SN2B':7,
            'SN2C':8}
    
    PROPNAME = {'APL':'APL',
                'CUR':r'\kappa',
                'THI':r'Z^{\prime}',
                'SN1A':r'S_{1A}',
                'SN1B':r'S_{1B}',
                'SN1C':r'S_{1C}',
                'SN2A':r'S_{2A}',
                'SN2B':r'S_{2B}',
                'SN2C':r'S_{2C}'}
    
    fig = plt.figure(figsize = (10, 10))
    
    gs = GS(5, 4, figure = fig, height_ratios = (1,1,1,1,0.25))
    
    ticks = np.linspace(10, 50, 5)
    
    n = 0
    for i in range(4):
        
        for j in range(4):
            
            data = PVAL[n]
            
            ax = fig.add_subplot(gs[i:i+1,j:j+1])
            
            if j == 0:
                
                ax.set_ylabel(r'$\mathbf{\nu}$ / MHz',
                              fontsize = 12, fontweight = 'bold')
                ax.set_yticks(ticks)
            
            elif j == 3:
                
                twin = ax.twinx()
                twin.set_ylabel(f'{MEMBRANES[i]}',
                                fontsize = 16, fontweight = 'bold',
                                rotation = 0, labelpad = 27)
                twin.set_yticks([])
                ax.set_yticks([])
                
            else:
                ax.set_yticks([])
                
            if i == 3:
                ax.set_xlabel('P / bar',
                              fontsize = 12, fontweight = 'bold')
                ax.set_xticks(ticks)
            
            elif i == 0:
                ax.set_title(f'{MEMBRANES[j]}',
                             fontsize = 16, fontweight = 'bold')
                ax.set_xticks([])
                
            else:
                ax.set_xticks([])
            
            
            img = ax.pcolormesh(PRES, FREQ, data[:,:,PROP[P.upper()]],
                                cmap = 'inferno_r',
                                vmin = MINPV[PROP[P.upper()]], vmax = MAXPV[PROP[P.upper()]])
            
            n += 1
    
    
    cbar_ticks = np.linspace(MINPV[PROP[P.upper()]], MAXPV[PROP[P.upper()]], 4) 
    ax_cbar    = fig.add_subplot(gs[-1,:])
    
    cbar = fig.colorbar(img, cax = ax_cbar,
                        ticks = cbar_ticks,
                        location = 'top')
    cbar.set_label(r'$\mathbf{p_{%s}}$' % (PROPNAME[P.upper()]),
                   fontsize = 16, fontweight = 'bold',
                   rotation = 0, labelpad = -50)
    
    cbar.ax.tick_params(top = False, labeltop = False,
                        bottom = True, labelbottom = True)
    cbar.ax.set_xticklabels(['{:4.2E}'.format(val) for val in cbar_ticks],
                            fontsize = 12)
    
    fig.tight_layout()
    
    plt.savefig(f'pval_{P.upper()}.svg',
                dpi = 600)

if __name__ == '__main__':
    
    DATA = []
    
    with open('./PVALS.pk', 'rb') as f:
        PVAL = pk.load(f)
    f.close()
    
    for membrane in MEMBRANES:
        
        with open(f'./{membrane}/PRES.pk', 'rb') as f:
            PRES = pk.load(f)
        f.close()
    
        with open(f'./{membrane}/FREQ.pk', 'rb') as f:
            FREQ = pk.load(f)
        f.close()
    
        with open(f'./{membrane}/DATA.pk', 'rb') as f:
            DATA.append(pk.load(f))
        f.close()
    
    DATA = np.array(DATA)
    
    DIFF = []
    size_data = np.loadtxt('size_data.txt')
    

    for i in range(4):
        
        data_i = DATA[i]
        
        for j in range(4):
            
            data_j = DATA[j]
            
            diff = data_j - data_i
            
            DIFF.append(diff)
            
    
    DIFF = np.array(DIFF)
    PVAL = PVAL.reshape(DIFF.shape)
    
    MEAN_DIFF = []
    STDV_DIFF = []
    
    dummy = []
    
    n = 0
    for i in range(4):
        for j in range(i+1, 4):
            
            dummy.append(n)
            n+=1
    
    dummy = np.array(dummy)
    
    MAX = []
    MIN = []
    
    for n in range(9):
        
        diffdata = DIFF[dummy,:,:,n,0]
        pvaldata = PVAL[dummy,:,:,n,0]
        
        idx_pf = np.logical_and(PRES < 35, FREQ > 15)[None,:,:]
        idx_pv = pvaldata < 1e20
        
        idx = np.logical_and(idx_pf, idx_pv)
        
        m = np.mean(diffdata[idx])
        s = np.std(diffdata[idx])
        
            
        MAX.append(10*s)
        MIN.append(-10*s)
    
    MEAN_PV = np.mean(PVAL, axis = (0,1,2))
    STDV_PV = np.std(PVAL,  axis = (0,1,2))
    
    MAXPV = [sig_level for m, s in zip(MEAN_PV, STDV_PV)]
    MINPV = [0 for m, s in zip(MEAN_PV, STDV_PV)]
    
    
    Graph('apl')
    Graph('cur')
    Graph('thi')
    Graph('sn1a')
    Graph('sn1b')
    Graph('sn1c')
    Graph('sn2a')
    Graph('sn2b')
    Graph('sn2c')
    
