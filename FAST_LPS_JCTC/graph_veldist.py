# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 13:53:55 2023

@author: alexandre
"""

import os
import numpy             as np
import pickle            as pk
import cblind            as cb
import statsmodels.api   as sm
import matplotlib.pyplot as plt

from scipy.integrate import simps

from statsmodels.tsa.stattools import acf

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


if __name__ == '__main__':
    
    for replica in range(2,4,1):
        
        for folder in ['NoFast', 'HalfFast', 'AllFast']:
            
            fullpath = f'Rep{replica}/{folder}'
        
            with open(f'{fullpath}/velocities.pk', 'rb') as f:
                VELS = pk.load(f)
            f.close()
        
            dt = 3000000/VELS.VELS.shape[0]
            W  = int(np.ceil(100000/dt))
        
            idx_0 = int(np.ceil(2000000/dt))
            idx_1 = int(np.ceil(3000000/dt))
        
            BOUND  = VELS.BOUND[idx_0:idx_1]
            BOUNDF = VELS.BOUNDF[idx_0:idx_1]
            VELS   = VELS.VELS[idx_0:idx_1]
        
            colors, linestyles = cb.Colorplots().cblind(2)
        
            fig, ax = plt.subplots(figsize = (7,7))
        
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
        
        
                ax.plot(kde_bound_remp.support, kde_bound_remp.density,
                        color = colors[0], label = 'Bound LPS')
                ax.plot(kde_unbound_remp.support, kde_unbound_remp.density,
                        color = colors[1], label = 'Unbound LPS')
        
        
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
        
                ax.plot(kde_bound_femp.support, kde_bound_femp.density,
                        color = colors[0], ls = '--',  label = 'Bound Fast LPS')
                ax.plot(kde_unbound_femp.support, kde_unbound_femp.density,
                        color = colors[1], ls = '--', label = 'Unbound Fast LPS')
        
        
            ax.set_xlabel(r'$\mathbf{v_{3D}}$ / $\mathbf{\AA \cdot ps^{-1}}$',
                          fontsize = 20, fontweight = 'bold')
        
            ax.set_ylabel(r'$\mathbf{P(v_{3D})}$',
                          fontsize = 20, fontweight = 'bold')
            ax.tick_params(labelsize = 18)
        
            ax.set_ylabel(r'$\mathbf{P(v_{3D})}$',
                          fontsize = 20, fontweight = 'bold')
            ax.tick_params(labelsize = 18)
        
            ax.set_xlim(-0.005, 1)
        
            ax.legend(fontsize = 16)
        
            fig.tight_layout()
        
            plt.savefig(f'{fullpath}/veldist.png')



