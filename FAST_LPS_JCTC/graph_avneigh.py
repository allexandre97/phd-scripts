# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 14:37:36 2023

@author: alexandre
"""

import numpy  as np
import pickle as pk
import matplotlib.pyplot as plt

from WHAM.lib                  import timeseries
from scipy.integrate           import simpson
from statsmodels.tsa.stattools import acf

def GenSubsamples(DATA, sepindex):
    
    out    = []
    checks = []
    
    idx0 = 0
    
    while len(checks) < len(DATA):
        
        subsample = []
        i = idx0
        while i < len(DATA):
            
            subsample.append(DATA[i])
            checks.append(i)
            i += sepindex
            
        out.append(subsample)
        idx0 += 1
        
    return out

def FDBins(X,
           Min = None, Max = None):
    
    q25, q75 = np.quantile(X, [0.25, 0.75])
    
    iqr = q75 - q25
    
    bw = 2*iqr/(X.shape[0]**(1/3))
    
    if Min == None:
        MIN = X.min()
    else:
        MIN = Min
    if Max == None:
        MAX = X.max()
    else:
        MAX = Max
    
    bins = np.arange(MIN, MAX + bw, bw)
    
    return bins

def ExpectedValue(X, PDF_x, dx):
    
    Series = [X[0] * PDF_x[0]]
    
    for n, (x, p) in enumerate(zip(X[1:], PDF_x[1:])):
        
        Series.append(x*p+ Series[n])
    
    return np.array(Series) * dx

def StandardDev(X, E, PDF_x, dx):
    
    EX2 = ExpectedValue(X**2, PDF_x, dx)
    
    return np.sqrt(EX2 - E**2), np.sqrt( (EX2[-1]-EX2[0]) - (E[-1]-E[0])**2 )


def MovAvFilter(data, win):
    
    windows = np.lib.stride_tricks.sliding_window_view(data, win)
    
    out = np.zeros((windows.shape[0], 2))
    
    out[:,0] = np.mean(windows, axis = 1)
    out[:,1] = np.std(windows, axis = 1)
    
    return out

def BuildViolinPlots(DATA, ax, i0 = 0, i1 = -1):
    
    DATA = np.array(DATA)[:,i0:i1]
    parts = ax.violinplot(DATA.T, showmeans=False, showmedians=False,
                             showextrema=False)
    
    for n, body in enumerate(parts['bodies']):
        
        body.set_color(colors[n])
        body.set_edgecolor('k')
        body.set_alpha(0.85)
    
    q1, ms, q3 = np.percentile(DATA, [25, 50, 75], axis =1)
    means = np.mean(DATA, axis = 1)
    w_min, w_max = means-1.5*(q3-q1), means+1.5*(q3-q1)

    inds = np.arange(1, len(ms) + 1)

    ax.scatter(inds, means, marker='X', color='gold', s=45, zorder=3)
    ax.vlines(inds, q1, q3, color='seagreen', linestyle='-', lw=5)
    ax.vlines(inds, w_min, w_max, color='seagreen', linestyle='-', lw=2)
    
def ACF(DATA, dataname, path):
    
    fig, ax = plt.subplots(figsize = (7,7))
    
    nlags = int(np.ceil(1000/(X[1]-X[0])))
    
    tau = np.linspace(0, 1000, nlags+1)
    
    G    = []
    SUBS = []
    for n in range(3):
        
        a, alpha = acf(DATA[n][i0:], nlags = nlags, alpha = 0.05)
        
        stat_ineff = timeseries.statisticalInefficiency(DATA[n][i0:])
        
        subsamples = GenSubsamples(DATA[n][i0:], int(stat_ineff))
        
        SUBS.append(subsamples)
        
        
        
        g = simpson(a, tau, dx = tau[1]-tau[0])
        G.append(g)
        
        ax.plot(tau, a, lw = 3,
                c = colors[n], alpha = 0.75,
                label = folderlabel[FOLDERS[n]])
        ax.fill_between(tau, alpha[:,0], alpha[:,1],
                        color = colors[n], edgecolor = None,
                        alpha = 0.25)
    
    ax.set_xlabel(r'$\mathbf{\tau}$ / ns', fontsize = 20, fontweight = 'bold')
    ax.set_ylabel(rf'ACF({dataname})', fontsize = 20, fontweight = 'bold')
    ax.tick_params(labelsize = 18)
    ax.legend(fontsize = 16)
    fig.tight_layout()
    
    plt.savefig(path,
                format = 'svg',
                dpi    = 600)
    
    return G, SUBS

def GraphSubsamples(Subs, DATA,
                    dataname,
                    Min, Max,
                    path,
                    bins = None):
    
    fig, ax = plt.subplots(3, 1,
                           sharex = True, sharey = True,
                           figsize = (7, 9))
    
    EXPVAL = []
    SIGMA  = []
    
    for n in range(3):
        
        
        ax[n].set_title(folderlabel[FOLDERS[n]],
                        fontsize = 20, fontweight = 'bold')
        
        ax2 = ax[n].twinx()
        
        subs = Subs[n]
        
        if type(bins) == type(None):
            BINS = FDBins(DATA[n][i0:], Min = Min, Max = Max)
        else:
            BINS = bins
        
        B = 0.5*(BINS[1:] + BINS[:-1])
        
        Ex, E, Sg, S = [], [], [], []
        
        for s in subs:
            
            N, _, _ = ax[n].hist(s, bins = BINS, density = True,
                                 color = 'royalblue', edgecolor = 'k', alpha = 0.05)
            
            expval = ExpectedValue(B, N, B[1] - B[0])
            sigma, s = StandardDev(B, expval, N, B[1] - B[0])
            
            Ex.append(expval)
            E.append(expval[-1] - expval[0])
            Sg.append(sigma)
            S.append(s)
            
        
        Ex = np.array(Ex)
        Sg = np.array(Sg)
        
        EXPVAL.append(E)
        SIGMA.append(S)
        
        Ex_m = np.mean(Ex, axis = 0)
        Ex_s = np.std(Ex, axis = 0)
        
        Sg_m = np.mean(Sg, axis = 0)
        Sg_s = np.std(Sg, axis = 0)
        
        
        ax2.plot(B, Ex_m,
                 c = 'gold', alpha = 1, label = r'$\mathbf{\mu (\mathbb{E}(%s))}$ = %.3f $\pm$ %.3f' % (dataname, np.mean(E), 3*np.std(E)))
        ax2.fill_between(B, 
                         Ex_m - Ex_s*3, Ex_m + Ex_s*3,
                         color = 'gold', edgecolor = None, alpha = 0.15)
        
        ax2.plot(B, Sg_m,
                 c = 'seagreen', ls = '--', alpha = 1, label = r'$\mathbf{\mu (\mathbb{\sigma}(%s))}$ = %.3f $\pm$ %.3f' % (dataname, np.mean(S), 3*np.std(S)))
        ax2.fill_between(B, 
                         Sg_m - Sg_s*3, Sg_m + Sg_s*3,
                         color = 'seagreen', edgecolor = None, alpha = 0.15)
        
        
        ax2.legend(fontsize = 12, loc = 'upper left')
        ax2.set_ylim(0, Max)
        ax2.tick_params(labelsize = 18)
        ax2.set_ylabel(r'$\mathbf{%s}$' % dataname, fontsize = 20, fontweight = 'bold')
        
        ax[n].set_ylabel(r'P($\mathbf{%s}$)' % dataname, fontsize = 20, fontweight = 'bold')
        ax[n].tick_params(labelsize = 18)
    
    ax[-1].set_xlabel(r'$\mathbf{%s}$qq' % dataname, fontsize = 20, fontweight = 'bold')
    
    fig.tight_layout()
    
    plt.savefig(path,
                format = 'svg',
                dpi = 600)
    
    return EXPVAL, SIGMA

def GraphViolin(E,
                dataname,
                path):
    
    fig, ax = plt.subplots(figsize = (7, 7))
    
    parts = ax.violinplot(E, showmeans=False, showmedians=False,
                             showextrema=False)
    
    for n, body in enumerate(parts['bodies']):
        
        body.set_color(colors[n])
        body.set_edgecolor('k')
        body.set_alpha(0.85)
    
    q1, ms, q3 = [np.percentile(E[n], 25) for n in range(3)], [np.percentile(E[n], 50) for n in range(3)], [np.percentile(E[n], 75) for n in range(3)]
    q1 = np.array(q1)
    ms = np.array(ms)
    q3 = np.array(q3)
    means = [np.mean(E[n]) for n in range(3)]
    means = np.array(means)
    w_min, w_max = means-1.5*(q3-q1), means+1.5*(q3-q1)

    inds = np.arange(1, len(ms) + 1)

    ax.scatter(inds, means, marker='X', color='gold', s=45, zorder=3)
    ax.vlines(inds, q1, q3, color='seagreen', linestyle='-', lw=5)
    ax.vlines(inds, w_min, w_max, color='seagreen', linestyle='-', lw=2)
    
    ax.set_xticks([1,2,3])
    ax.set_xticklabels([folderlabel[FOLDERS[n]] for n in range(3)])
    ax.tick_params(labelsize = 18)
    ax.set_xlabel('System', fontsize = 20, fontweight = 'bold')
    ax.set_ylabel(r'$\mathbf{\mathbb{E}(%s)}$' % (dataname), fontsize = 20, fontweight = 'bold')

    fig.tight_layout()
    
    plt.savefig(path,
                format = 'svg',
                dpi = 600)
    
FOLDERS = ['NoFast',
           'HalfFast',
           'AllFast']

folderlabel = {'NoFast':'0% FastLPS',
               'HalfFast':'50% FastLPS',
               'AllFast':'100% FastLPS'}

colors = ['#648fff', '#fe6100', '#000000']


if __name__ == '__main__':
    
    FULLDATA_DENS = [[],[],[]]
    FULLDATA_NCLU = [[],[],[]]
    
    for rep in ['Rep1', 'Rep2', 'Rep3']:
    
        fig, ax = plt.subplots(1, 2,
                               sharey  = True,
                               figsize = (12, 7))
        
        fig2, ax2 = plt.subplots(1, 2,
                                 sharey  = True,
                                 figsize = (12, 7))
        
        DATA  = []
        DATA2 = []
        
        for n, fast in enumerate(['NoFast', 'HalfFast', 'AllFast']):
            
            with open(f'{rep}/{fast}/clusters_avnei_0ns_5000ns_2.pk', 'rb') as f:
                data = pk.load(f)
            f.close()
            
            with open(f'{rep}/{fast}/clusters_nclust_0ns_5000ns.pk', 'rb') as f:
                data2 = pk.load(f)
            f.close()
            
            X = data[:,0]/1000
            Y = data[:,1]
            
            X2 = data2[:,0]/1000
            Y2 = data2[:,1]
            
            DATA.append(Y)
            DATA2.append(Y2)
            
            win = int(np.ceil(10/(X[1]-X[0])))
            
            Y_filtered  = MovAvFilter(Y,  win)
            Y2_filtered = MovAvFilter(Y2, win)
            
            if win % 2 == 0:
                win += 1
            
            ax[0].plot(X[:Y_filtered.shape[0]], Y_filtered[:,0],
                       c = colors[n], lw = 3, alpha = 0.75,
                       label = folderlabel[fast])
            ax[0].fill_between(X[:Y_filtered.shape[0]],
                               Y_filtered[:,0] + Y_filtered[:,1]*3,
                               Y_filtered[:,0] - Y_filtered[:,1]*3,
                               color = colors[n], alpha = 0.25)
            
            ax2[0].plot(X2[:Y2_filtered.shape[0]], Y2_filtered[:,0],
                        c = colors[n], lw = 3, alpha = 0.75,
                        label = folderlabel[fast])
            ax2[0].fill_between(X2[:Y2_filtered.shape[0]],
                                Y2_filtered[:,0] + Y2_filtered[:,1]*3,
                                Y2_filtered[:,0] - Y2_filtered[:,1]*3,
                                color = colors[n], alpha = 0.25)
        
        
        i0 = int(np.ceil(2000/(X[1]-X[0])))
        
        for _ in range(3):
            
            for tmp1 in range(len(DATA[_])):
                FULLDATA_DENS[_].append(DATA[_][tmp1])
            for tmp2 in range(len(DATA2[_])):
                FULLDATA_NCLU[_].append(DATA2[_][tmp2])
        
# =============================================================================
#         BuildViolinPlots(DATA,  ax[1], i0)
#         BuildViolinPlots(DATA2, ax2[1], i0)
#         
#         ax[0].set_xlabel('Time / ns', fontsize = 20, fontweight = 'bold')
#         ax[0].set_ylabel(r'$\mathbf{\rho}$', fontsize = 20, fontweight = 'bold')
#         ax[0].tick_params(labelsize = 18)
#         ax[0].legend(fontsize = 16)
#         
#         ax[1].set_xlabel('System', fontsize = 20, fontweight = 'bold')
#         ax[1].set_xticks([1, 2, 3])
#         ax[1].set_xticklabels(list(folderlabel.values()))
#         ax[1].tick_params(labelsize = 18)
#         
#         fig.tight_layout()
#         
#         plt.savefig(f'{rep}/cluster_density.svg',
#                     format = 'svg',
#                     dpi    = 600)
#         
#         
#         ax2[0].set_xlabel('Time / ns', fontsize = 20, fontweight = 'bold')
#         ax2[0].set_ylabel('Number of Clusters', fontsize = 20, fontweight = 'bold')
#         ax2[0].tick_params(labelsize = 18)
#         ax2[0].legend(fontsize = 16)
#         
#         ax2[1].set_xlabel('System', fontsize = 20, fontweight = 'bold')
#         ax2[1].set_xticks([1, 2, 3])
#         ax2[1].set_xticklabels(list(folderlabel.values()))
#         ax2[1].tick_params(labelsize = 18)
#         
#         fig2.tight_layout()
#         
#         plt.savefig(f'{rep}/n_clusters.svg',
#                     format = 'svg',
#                     dpi    = 600)
#         
#         
#         G_avneig, Subs_avneig = ACF(DATA,  r'$\mathbf{\rho}$', f'{rep}/acf_density.svg')
#         G_nclust, Subs_nclust = ACF(DATA2, 'n. clusters', f'{rep}/acf_nclusters.svg')
#         
#         
#         E_AVN, S_AVN = GraphSubsamples(Subs_avneig, DATA,  r'\rho', 0.4, 1.0, 
#                                        f'{rep}/bootstrap_density.svg')
#         E_NCL, S_NCL = GraphSubsamples(Subs_nclust, DATA2, 'n. clusters', 5.0, 30,
#                                        f'{rep}/bootstrap_nclusters.svg',
#                                        np.arange(5, 30, 1))
#         
#         
#         GraphViolin(E_AVN, r'\rho', f'{rep}/violin_density.svg')
#         GraphViolin(E_NCL, 'n. clusters', f'{rep}/violin_nclusters.svg')
# =============================================================================
    
    def bootstrap_array(data, num_subsamples, subsample_size=None):
        """
        Bootstrap the given array by generating subsamples with replacement.
    
        Parameters:
            data (numpy array): The original array to be bootstrapped.
            num_subsamples (int): Number of subsamples to generate.
            subsample_size (int or None): Size of each subsample. If None, use the size of the original data.
    
        Returns:
            list: A list of numpy arrays, each representing a subsample.
        """
        subsamples = []
        n = len(data)
    
        if subsample_size is None:
            subsample_size = n  # Default to the size of the original data
    
        for _ in range(num_subsamples):
            # Generate a random sample with replacement
            subsample = np.random.choice(data, subsample_size, replace=True)
            subsamples.append(subsample)
    
        return subsamples
    
    AVG_DNS = []
    AVG_NCL = []
    
    for n in range(3):
        
        boots_dns = np.mean(bootstrap_array(FULLDATA_DENS[n], 1000), axis = 0)
        boots_ncl = np.mean(bootstrap_array(FULLDATA_NCLU[n], 1000), axis = 0)
        
        AVG_DNS.append(boots_dns)
        AVG_NCL.append(boots_ncl)
    
    GraphViolin(AVG_DNS, r'\rho', './avg_violin_density.svg')
    GraphViolin(AVG_NCL, 'n. clusters', './violin_nclusters.svg')

    
    
    
    
    
    
    
        