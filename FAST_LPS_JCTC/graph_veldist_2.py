# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 13:53:55 2023

@author: alexandre
"""

import numpy             as np
import pickle            as pk
import statsmodels.api   as sm
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d


class Velocities:

    def __init__(self,
                 BOUND,
                 BOUNDF,
                 VELS):

        self.BOUND = BOUND
        self.BOUNDF = BOUNDF
        self.VELS  = VELS

def GetBoundUnbound(BOUND, VELS):

    wherebound   = np.where(BOUND[:-1])
    whereunbound = np.where(np.logical_not(BOUND[:-1]))

    vels_bound   = VELS[wherebound[0], wherebound[1]]
    vels_unbound = VELS[whereunbound[0], whereunbound[1]]

    v3d_bound   = np.sqrt(vels_bound[:,0]**2 + vels_bound[:,1]**2 + vels_bound[:,2]**2)
    v3d_unbound = np.sqrt(vels_unbound[:,0]**2 + vels_unbound[:,1]**2 + vels_unbound[:,2]**2)

    return v3d_bound, v3d_unbound


folderlabel = {'NoFast':'0% FastLPS',
               'HalfFast':'50% FastLPS',
               'AllFast':'100% FastLPS'}

if __name__ == '__main__':
    
    
    colors = ['#648fff', '#fe6100', '#000000', '#dc267f']
    
    FULLDATA = []
    
    for replica in range(1, 4, 1):
        
        REP = []
        
        fig, ax = plt.subplots(1, 2, figsize = (12,4),
                               sharex = True, sharey = True)
        
        for n, folder in enumerate(['NoFast', 'HalfFast', 'AllFast']):
            
            REMP, FEMP = [], []
            
            with open(f'Rep{replica}/{folder}/velocities.pk', 'rb') as f:
                VELS = pk.load(f)
            f.close()
        
            dt = 5000000/VELS.VELS.shape[0]
            W  = int(np.ceil(100000/dt))
        
            idx_0 = int(np.ceil(4000000/dt))
            idx_1 = int(np.ceil(5000000/dt))
        
            BOUND  = VELS.BOUND[idx_0:idx_1]
            BOUNDF = VELS.BOUNDF[idx_0:idx_1]
            VELS   = VELS.VELS[idx_0:idx_1]
        
        
            if (VELS[:,0,:,:] != 0).any():
        
                v3d_bound_remp, v3d_unbound_remp = GetBoundUnbound(BOUND,  VELS[:,0,:,:])
        
                with open('v3d_bound_remp.pk', 'wb') as f:
                    pk.dump(v3d_bound_remp, f, protocol = 4)
                f.close()
        
                with open('v3d_unbound_remp.pk', 'wb') as f:
                    pk.dump(v3d_unbound_remp, f, protocol = 4)
                f.close()
        
                kde_bound_remp = sm.nonparametric.KDEUnivariate(v3d_bound_remp)
                kde_bound_remp.fit(gridsize = 10000)
        
                kde_unbound_remp = sm.nonparametric.KDEUnivariate(v3d_unbound_remp)
                kde_unbound_remp.fit(gridsize = 10000)
                
                REMP.append([[kde_bound_remp.support, kde_bound_remp.density],
                             [kde_unbound_remp.support, kde_unbound_remp.density]])
        
                ax[0].plot(kde_bound_remp.support, kde_bound_remp.density,
                           color = colors[n], ls = '-', label = f'{folderlabel[folder]}, Bound LPS')
                ax[0].plot(kde_unbound_remp.support, kde_unbound_remp.density,
                           color = colors[n], ls = '--', label = f'{folderlabel[folder]}, Unbound LPS')
        
        
            if (VELS[:,1,:,:] != 0).any():
        
                v3d_bound_femp, v3d_unbound_femp = GetBoundUnbound(BOUNDF, VELS[:,1,:,:])
        
                with open('v3d_bound_femp.pk', 'wb') as f:
                    pk.dump(v3d_bound_femp, f, protocol = 4)
                f.close()
        
                with open('v3d_unbound_femp.pk', 'wb') as f:
                    pk.dump(v3d_unbound_femp, f, protocol = 4)
                f.close()
        
                kde_bound_femp = sm.nonparametric.KDEUnivariate(v3d_bound_femp)
                kde_bound_femp.fit(gridsize = 10000)
        
                kde_unbound_femp = sm.nonparametric.KDEUnivariate(v3d_unbound_femp)
                kde_unbound_femp.fit(gridsize = 10000)
                
                FEMP.append([[kde_bound_femp.support, kde_bound_femp.density],
                             [kde_unbound_femp.support, kde_unbound_femp.density]])
        
                ax[1].plot(kde_bound_femp.support, kde_bound_femp.density,
                        color = colors[n], ls = '-', label = f'{folderlabel[folder]}, Bound Fast LPS')
                ax[1].plot(kde_unbound_femp.support, kde_unbound_femp.density,
                        color = colors[n], ls = '--', label = f'{folderlabel[folder]}, Unbound Fast LPS')
            
            REP.append([REMP, FEMP])
        
        FULLDATA.append(REP)
        
        
        ax[0].set_xlabel(r'$\mathbf{v_{3D}}$ / $\mathbf{\AA \cdot ps^{-1}}$',
                      fontsize = 20, fontweight = 'bold')
    
        ax[0].set_ylabel(r'$\mathbf{P(v_{3D})}$',
                      fontsize = 20, fontweight = 'bold')
        ax[0].tick_params(labelsize = 18)
    
        ax[0].legend(fontsize = 14)
    
        ax[1].set_xlabel(r'$\mathbf{v_{3D}}$ / $\mathbf{\AA \cdot ps^{-1}}$',
                      fontsize = 20, fontweight = 'bold')
    
        ax[1].tick_params(labelsize = 18)
    
        ax[1].set_xlim(-0.005, 1)
        ax[1].set_ylim(-0.05, 15)
    
        ax[1].legend(fontsize = 14)
    
        fig.tight_layout()
        
        plt.savefig(f'./Rep{replica}/veldist_final.svg',
                    format = 'svg',
                    dpi = 600)
        
    

def InterpolateFulldata(FULLDATA):
    
    colors = ['#648fff', '#fe6100', '#000000', '#dc267f']
    
    X = np.unique(np.concatenate((#REMP UNBOUND NOFAST
                                  FULLDATA[0][0][0][0][0][0],
                                  FULLDATA[1][0][0][0][0][0],
                                  FULLDATA[2][0][0][0][0][0],
                                  #REMP UNBOUND HALFFAST
                                  FULLDATA[0][1][0][0][0][0],
                                  FULLDATA[1][1][0][0][0][0],
                                  FULLDATA[2][1][0][0][0][0],
                                  #REMP BOUND NOFAST
                                  FULLDATA[0][0][0][0][1][0],
                                  FULLDATA[1][0][0][0][1][0],
                                  FULLDATA[2][0][0][0][1][0],
                                  #REMP BOUND HALFFAST
                                  FULLDATA[0][1][0][0][1][0],
                                  FULLDATA[1][1][0][0][1][0],
                                  FULLDATA[2][1][0][0][1][0],
                                  
                                  #FAST REMP UNBOUND HALFFAST
                                  FULLDATA[0][1][1][0][0][0],
                                  FULLDATA[1][1][1][0][0][0],
                                  FULLDATA[2][1][1][0][0][0],
                                  #FAST REMP UNBOUND ALLFAST
                                  FULLDATA[0][2][1][0][0][0],
                                  FULLDATA[1][2][1][0][0][0],
                                  FULLDATA[2][2][1][0][0][0],
                                  #FAST REMP BOUND HALFFAST
                                  FULLDATA[0][1][1][0][1][0],
                                  FULLDATA[1][1][1][0][1][0],
                                  FULLDATA[2][1][1][0][1][0],
                                  #FAST REMP BOUND ALLFAST
                                  FULLDATA[0][2][1][0][1][0],
                                  FULLDATA[1][2][1][0][1][0],
                                  FULLDATA[2][2][1][0][1][0])))
    
    Y1_REMP_UNBOUND_NOFAST = interp1d(FULLDATA[0][0][0][0][0][0],
                                      FULLDATA[0][0][0][0][0][1],
                                      kind = 'linear', fill_value = 'extrapolate')(X)
    Y2_REMP_UNBOUND_NOFAST = interp1d(FULLDATA[1][0][0][0][0][0],
                                      FULLDATA[1][0][0][0][0][1],
                                      kind = 'linear', fill_value = 'extrapolate')(X)
    Y3_REMP_UNBOUND_NOFAST = interp1d(FULLDATA[2][0][0][0][0][0],
                                      FULLDATA[2][0][0][0][0][1],
                                      kind = 'linear', fill_value = 'extrapolate')(X)
    
    Y_REMP_UNBOUND_NOFAST = np.array([np.mean([Y1_REMP_UNBOUND_NOFAST,
                                               Y2_REMP_UNBOUND_NOFAST,
                                               Y3_REMP_UNBOUND_NOFAST], axis = 0),
                                      np.std([Y1_REMP_UNBOUND_NOFAST,
                                              Y2_REMP_UNBOUND_NOFAST,
                                              Y3_REMP_UNBOUND_NOFAST], axis = 0)])
    
    Y1_REMP_BOUND_NOFAST = interp1d(FULLDATA[0][0][0][0][1][0],
                                    FULLDATA[0][0][0][0][1][1],
                                    kind = 'linear', fill_value = 'extrapolate')(X)
    Y2_REMP_BOUND_NOFAST = interp1d(FULLDATA[1][0][0][0][1][0],
                                    FULLDATA[1][0][0][0][1][1],
                                    kind = 'linear', fill_value = 'extrapolate')(X)
    Y3_REMP_BOUND_NOFAST = interp1d(FULLDATA[2][0][0][0][1][0],
                                    FULLDATA[2][0][0][0][1][1],
                                    kind = 'linear', fill_value = 'extrapolate')(X)
    
    Y_REMP_BOUND_NOFAST = np.array([np.mean([Y1_REMP_BOUND_NOFAST,
                                             Y2_REMP_BOUND_NOFAST,
                                             Y3_REMP_BOUND_NOFAST], axis = 0),
                                    np.std([Y1_REMP_BOUND_NOFAST,
                                            Y2_REMP_BOUND_NOFAST,
                                            Y3_REMP_BOUND_NOFAST], axis = 0)])
    

    Y1_REMP_UNBOUND_HALFFAST = interp1d(FULLDATA[0][1][0][0][0][0],
                                        FULLDATA[0][1][0][0][0][1],
                                        kind = 'linear', fill_value = 'extrapolate')(X)
    Y2_REMP_UNBOUND_HALFFAST = interp1d(FULLDATA[1][1][0][0][0][0],
                                        FULLDATA[1][1][0][0][0][1],
                                        kind = 'linear', fill_value = 'extrapolate')(X)
    Y3_REMP_UNBOUND_HALFFAST = interp1d(FULLDATA[2][1][0][0][0][0],
                                        FULLDATA[2][1][0][0][0][1],
                                        kind = 'linear', fill_value = 'extrapolate')(X)
    
    Y_REMP_UNBOUND_HALFFAST = np.array([np.mean([Y1_REMP_UNBOUND_HALFFAST,
                                                 Y2_REMP_UNBOUND_HALFFAST,
                                                 Y3_REMP_UNBOUND_HALFFAST], axis = 0),
                                        np.std([Y1_REMP_UNBOUND_HALFFAST,
                                                Y2_REMP_UNBOUND_HALFFAST,
                                                Y3_REMP_UNBOUND_HALFFAST], axis = 0)])
    
    Y1_REMP_BOUND_HALFFAST = interp1d(FULLDATA[0][1][0][0][1][0],
                                      FULLDATA[0][1][0][0][1][1],
                                      kind = 'linear', fill_value = 'extrapolate')(X)
    Y2_REMP_BOUND_HALFFAST = interp1d(FULLDATA[1][1][0][0][1][0],
                                      FULLDATA[1][1][0][0][1][1],
                                      kind = 'linear', fill_value = 'extrapolate')(X)
    Y3_REMP_BOUND_HALFFAST = interp1d(FULLDATA[2][1][0][0][1][0],
                                      FULLDATA[2][1][0][0][1][1],
                                      kind = 'linear', fill_value = 'extrapolate')(X)
    
    Y_REMP_BOUND_HALFFAST = np.array([np.mean([Y1_REMP_BOUND_HALFFAST,
                                               Y2_REMP_BOUND_HALFFAST,
                                               Y3_REMP_BOUND_HALFFAST], axis = 0),
                                      np.std([Y1_REMP_BOUND_HALFFAST,
                                              Y2_REMP_BOUND_HALFFAST,
                                              Y3_REMP_BOUND_HALFFAST], axis = 0)])
    
        

    Y1_FEMP_UNBOUND_HALFFAST = interp1d(FULLDATA[0][1][1][0][0][0],
                                        FULLDATA[0][1][1][0][0][1],
                                        kind = 'linear', fill_value = 'extrapolate')(X)
    Y2_FEMP_UNBOUND_HALFFAST = interp1d(FULLDATA[1][1][1][0][0][0],
                                        FULLDATA[1][1][1][0][0][1],
                                        kind = 'linear', fill_value = 'extrapolate')(X)
    Y3_FEMP_UNBOUND_HALFFAST = interp1d(FULLDATA[2][1][1][0][0][0],
                                        FULLDATA[2][1][1][0][0][1],
                                        kind = 'linear', fill_value = 'extrapolate')(X)
    
    Y_FEMP_UNBOUND_HALFFAST = np.array([np.mean([Y1_FEMP_UNBOUND_HALFFAST,
                                                 Y2_FEMP_UNBOUND_HALFFAST,
                                                 Y3_FEMP_UNBOUND_HALFFAST], axis = 0),
                                        np.std([Y1_FEMP_UNBOUND_HALFFAST,
                                                Y2_FEMP_UNBOUND_HALFFAST,
                                                Y3_FEMP_UNBOUND_HALFFAST], axis = 0)])
    
    Y1_FEMP_BOUND_HALFFAST = interp1d(FULLDATA[0][1][1][0][1][0],
                                      FULLDATA[0][1][1][0][1][1],
                                      kind = 'linear', fill_value = 'extrapolate')(X)
    Y2_FEMP_BOUND_HALFFAST = interp1d(FULLDATA[1][1][1][0][1][0],
                                      FULLDATA[1][1][1][0][1][1],
                                      kind = 'linear', fill_value = 'extrapolate')(X)
    Y3_FEMP_BOUND_HALFFAST = interp1d(FULLDATA[2][1][1][0][1][0],
                                      FULLDATA[2][1][1][0][1][1],
                                      kind = 'linear', fill_value = 'extrapolate')(X)
    
    Y_FEMP_BOUND_HALFFAST = np.array([np.mean([Y1_FEMP_BOUND_HALFFAST,
                                               Y2_FEMP_BOUND_HALFFAST,
                                               Y3_FEMP_BOUND_HALFFAST], axis = 0),
                                      np.std([Y1_FEMP_BOUND_HALFFAST,
                                              Y2_FEMP_BOUND_HALFFAST,
                                              Y3_FEMP_BOUND_HALFFAST], axis = 0)])
    
    Y1_FEMP_UNBOUND_ALLFAST = interp1d(FULLDATA[0][2][1][0][0][0],
                                       FULLDATA[0][2][1][0][0][1],
                                      kind = 'linear', fill_value = 'extrapolate')(X)
    Y2_FEMP_UNBOUND_ALLFAST = interp1d(FULLDATA[1][2][1][0][0][0],
                                       FULLDATA[1][2][1][0][0][1],
                                       kind = 'linear', fill_value = 'extrapolate')(X)
    Y3_FEMP_UNBOUND_ALLFAST = interp1d(FULLDATA[2][2][1][0][0][0],
                                       FULLDATA[2][2][1][0][0][1],
                                       kind = 'linear', fill_value = 'extrapolate')(X)
    
    Y_FEMP_UNBOUND_ALLFAST = np.array([np.mean([Y1_FEMP_UNBOUND_ALLFAST,
                                                Y2_FEMP_UNBOUND_ALLFAST,
                                                Y3_FEMP_UNBOUND_ALLFAST], axis = 0),
                                       np.std([Y1_FEMP_UNBOUND_ALLFAST,
                                               Y2_FEMP_UNBOUND_ALLFAST,
                                               Y3_FEMP_UNBOUND_ALLFAST], axis = 0)])
    
    Y1_FEMP_BOUND_ALLFAST = interp1d(FULLDATA[0][2][1][0][1][0],
                                     FULLDATA[0][2][1][0][1][1],
                                     kind = 'linear', fill_value = 'extrapolate')(X)
    Y2_FEMP_BOUND_ALLFAST = interp1d(FULLDATA[1][2][1][0][1][0],
                                     FULLDATA[1][2][1][0][1][1],
                                     kind = 'linear', fill_value = 'extrapolate')(X)
    Y3_FEMP_BOUND_ALLFAST = interp1d(FULLDATA[2][2][1][0][1][0],
                                     FULLDATA[2][2][1][0][1][1],
                                     kind = 'linear', fill_value = 'extrapolate')(X)
    
    Y_FEMP_BOUND_ALLFAST = np.array([np.mean([Y1_FEMP_BOUND_ALLFAST,
                                              Y2_FEMP_BOUND_ALLFAST,
                                              Y3_FEMP_BOUND_ALLFAST], axis = 0),
                                     np.std([Y1_FEMP_BOUND_ALLFAST,
                                             Y2_FEMP_BOUND_ALLFAST,
                                             Y3_FEMP_BOUND_ALLFAST], axis = 0)])
    
    fig, ax = plt.subplots(1, 2, figsize = (12,4),
                           sharex = True, sharey = True)
    
    
    ax[0].plot(X, Y_REMP_BOUND_NOFAST[0],
               c = colors[0], label = '0% FastLPS, Bound LPS')
    ax[0].plot(X, Y_REMP_UNBOUND_NOFAST[0],
               c = colors[0], ls = '--', label = '0% FastLPS, Unbound LPS')
    ax[0].fill_between(X, 
                       Y_REMP_BOUND_NOFAST[0] - Y_REMP_BOUND_NOFAST[1],
                       Y_REMP_BOUND_NOFAST[0] + Y_REMP_BOUND_NOFAST[1],
                       color = colors[0], edgecolor = None, alpha = 0.25)
    ax[0].fill_between(X, 
                       Y_REMP_UNBOUND_NOFAST[0] - Y_REMP_UNBOUND_NOFAST[1],
                       Y_REMP_UNBOUND_NOFAST[0] + Y_REMP_UNBOUND_NOFAST[1],
                       color = colors[0], edgecolor = None, alpha = 0.25)
    
    ax[0].plot(X, Y_REMP_BOUND_HALFFAST[0],
               c = colors[1], label = '50% FastLPS, Bound LPS')
    ax[0].plot(X, Y_REMP_UNBOUND_HALFFAST[0],
               c = colors[1], ls = '--', label = '50% FastLPS, Unbound LPS')
    ax[0].fill_between(X, 
                       Y_REMP_BOUND_HALFFAST[0] - Y_REMP_BOUND_HALFFAST[1],
                       Y_REMP_BOUND_HALFFAST[0] + Y_REMP_BOUND_HALFFAST[1],
                       color = colors[1], edgecolor = None, alpha = 0.25)
    ax[0].fill_between(X, 
                       Y_REMP_UNBOUND_HALFFAST[0] - Y_REMP_UNBOUND_HALFFAST[1],
                       Y_REMP_UNBOUND_HALFFAST[0] + Y_REMP_UNBOUND_HALFFAST[1],
                       color = colors[1], edgecolor = None, alpha = 0.25)


    ax[1].plot(X, Y_FEMP_BOUND_HALFFAST[0],
               c = colors[1], label = '50% FastLPS, Bound Fast LPS')
    ax[1].plot(X, Y_FEMP_UNBOUND_HALFFAST[0],
               c = colors[1], ls = '--', label = '50% FastLPS, Unbound Fast LPS')
    ax[1].fill_between(X, 
                       Y_FEMP_BOUND_HALFFAST[0] - Y_FEMP_BOUND_HALFFAST[1],
                       Y_FEMP_BOUND_HALFFAST[0] + Y_FEMP_BOUND_HALFFAST[1],
                       color = colors[1], edgecolor = None, alpha = 0.25)
    ax[1].fill_between(X, 
                       Y_FEMP_UNBOUND_HALFFAST[0] - Y_FEMP_UNBOUND_HALFFAST[1],
                       Y_FEMP_UNBOUND_HALFFAST[0] + Y_FEMP_UNBOUND_HALFFAST[1],
                       color = colors[1], edgecolor = None, alpha = 0.25)

    
    ax[1].plot(X, Y_FEMP_BOUND_ALLFAST[0],
               c = colors[2], label = '100% FastLPS, Bound Fast LPS')
    ax[1].plot(X, Y_FEMP_UNBOUND_ALLFAST[0],
               c = colors[2], ls = '--', label = '100% FastLPS, Unbound Fast LPS')
    ax[1].fill_between(X, 
                       Y_FEMP_BOUND_ALLFAST[0] - Y_FEMP_BOUND_ALLFAST[1],
                       Y_FEMP_BOUND_ALLFAST[0] + Y_FEMP_BOUND_ALLFAST[1],
                       color = colors[2], edgecolor = None, alpha = 0.25)
    ax[1].fill_between(X, 
                       Y_FEMP_UNBOUND_ALLFAST[0] - Y_FEMP_UNBOUND_ALLFAST[1],
                       Y_FEMP_UNBOUND_ALLFAST[0] + Y_FEMP_UNBOUND_ALLFAST[1],
                       color = colors[2], edgecolor = None, alpha = 0.25)    
    
    
    ax[0].set_xlabel(r'$\mathbf{v_{3D}}$ / $\mathbf{\AA \cdot ps^{-1}}$',
                  fontsize = 20, fontweight = 'bold')

    ax[0].set_ylabel(r'$\mathbf{P(v_{3D})}$',
                  fontsize = 20, fontweight = 'bold')
    ax[0].tick_params(labelsize = 18)

    ax[0].legend(fontsize = 14)

    ax[1].set_xlabel(r'$\mathbf{v_{3D}}$ / $\mathbf{\AA \cdot ps^{-1}}$',
                  fontsize = 20, fontweight = 'bold')

    ax[1].tick_params(labelsize = 18)

    ax[1].set_xlim(-0.005, 1)
    ax[1].set_ylim(-0.05, 15)

    ax[1].legend(fontsize = 14)

    fig.tight_layout()
    
    plt.savefig('avg_veldist_final.svg',
                format = 'svg',
                dpi = 600)

    
InterpolateFulldata(FULLDATA)