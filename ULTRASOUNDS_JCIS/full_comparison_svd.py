#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 17:13:15 2022

@author: alexandre
"""

import pickle
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import matplotlib.pyplot as plt



def Clusterize(dataframe, max_clusters):
    
    from kmodes.kprototypes import KPrototypes as kp
    from kneed import KneeLocator as KL
    
    comp = dataframe.pop('Composition')
    
    Cost = []
    for nc in range(2, max_clusters):
        
        try:
            kproto = kp(n_clusters = nc, init = 'Cao', verbose = 0, n_jobs = 10)
            clusters = kproto.fit_predict(dataframe, categorical=[4, 5, 6, 7, 12])
            cost = kproto.cost_
            Cost.append(cost)
        except:
            kproto = kp(n_clusters = nc, init = 'random', verbose = 0, n_jobs = 10)
            clusters = kproto.fit_predict(dataframe, categorical=[4, 5, 6, 7, 12])
            cost = kproto.cost_
            Cost.append(cost)
    
    kl = KL(range(2, max_clusters), Cost, curve = 'convex',  direction = 'decreasing')
    
    kproto = kp(n_clusters = kl.knee, init = 'Cao', n_jobs = 11)
    clusters = kproto.fit_predict(dataframe, categorical = [4, 5, 6, 7, 12])
    
    dataframe['Composition'] = comp
    dataframe['Cluster'] = clusters
    
    return dataframe, kl.knee
    
def Distributions(clusters, Property, colors, N):
    
    kde = []
    for dataframe in clusters:
        
        try:        
            kd = sm.nonparametric.KDEUnivariate(dataframe[Property])
            kd.fit()
            kde.append(kd)
        except:
            pass
            
    fig, ax = plt.subplots(figsize = (9,9))
    for i, kd in enumerate(kde):
        try:
            ax.plot(kd.support, kd.density*(clusters[i].shape[0]/N), label = 'Cluster '+str(i+1), color = colors[i])
            ax.fill_between(kd.support, kd.density*(clusters[i].shape[0]/N), alpha = 0.5, color = colors[i])
        except:
            pass
    ax.set_xlabel(Property, fontsize = 20, fontweight = 'bold')
    ax.set_ylabel('Probability Density', fontsize = 20, fontweight = 'bold')
    ax.tick_params(labelsize = 18)
    ax.legend(fontsize = 18)
    fig.tight_layout()
    
def Graphs(dataframe, nclust, N, total = True,):
        
    indices = []
    for index in range(nclust):
        indices.append(np.where(dataframe['Cluster'] == index)[0])
    
    clusters = []
    for index in indices:
        clusters.append(dataframe.iloc[index])
    
    
    nCP = []
    nFF = []

    for dfn in clusters:
        kk = []
        for cp in CP:
            try:
                n = dfn['Composition'].value_counts()[cp]
            except:
                n = 0
            kk.append(n)
        nCP.append(kk)
        kk= []
        for ff in ION:
            try:
                n = dfn['Ions'].value_counts()[ff]
            except:
                n = 0
            kk.append(n)
        nFF.append(kk)
    
    nCP = np.array(nCP)
    nFF = np.array(nFF)
    
    if total == True:
        fig, ax = plt.subplots(figsize = (12,9))
        ax.bar(range(1, nFF.shape[0]+1), nFF[:,0], color = 'royalblue', edgecolor = 'k', label = ION[0], alpha = 0.8)
        ax.bar(range(1, nFF.shape[0]+1), nFF[:,1], color = 'firebrick', edgecolor = 'k', label = ION[1], bottom = nFF[:,0], alpha = 0.8)
        ax.set_xlabel('Cluster Number', fontsize = 20, fontweight = 'bold')
        ax.set_ylabel('Number of Occurrences', fontsize = 20, fontweight = 'bold')        
        ax.set_xticks(range(1,nCP.shape[0]+1))
        ax.tick_params(labelsize = 18)
        ax.legend(fontsize = 18)
        fig.tight_layout()
        
    colors = ['royalblue', 'seagreen', 'firebrick', 'gold', 'indigo', 'turquoise', 'orangered', 'darkorange', 'navy', 'lightgreen']
        
    fig, ax = plt.subplots(figsize = (12, 9))
    
    ax.bar(range(1,nclust+1), nCP[:,0], label = CP[0], edgecolor = 'k', color = colors[0], alpha = 0.8)
    
    bottom = nCP[:,0].copy()
    
    for cp in range(1, len(CP)):
        ax.bar(range(1,nclust+1), nCP[:,cp], bottom = bottom, label = CP[cp], edgecolor = 'k', color = colors[cp], alpha = 0.8)
        bottom += nCP[:,cp].copy()
    
    ax.set_xlabel('Cluster Number', fontsize = 20, fontweight = 'bold')
    ax.set_ylabel('Number of Occurrences', fontsize = 20, fontweight = 'bold')
    
    ax.set_xticks(range(1,nCP.shape[0]+1))
    ax.tick_params(labelsize = 18)
    
    ax.legend(fontsize = 18)
    
    fig.tight_layout()
    
    for Property in dataframe.columns[0:4]:
        
        Distributions(clusters, Property, colors, N)
    
    NN = []
    
    for lipid in dataframe.columns[4:8]:
        N = []
        for cluster in clusters:
            
            n = len(np.where(cluster[lipid] != np.min(cluster[lipid]))[0])
            N.append(n)
        NN.append(N)
        
    NN = np.array(NN)
    
    fig, ax = plt.subplots(figsize = (12, 9))
    
    ax.bar(range(1, nclust+1), NN[0,:], label = CP[0], color = colors[0], edgecolor = 'k', alpha = 0.8)
    
    bottom = NN[0,:].copy()
    
    for i, lipid in enumerate(NN[1:]):
        ax.bar(range(1,nclust+1), lipid, bottom = bottom, label = CP[i+1], color = colors[i+1], edgecolor = 'k', alpha = 0.8)
        bottom += lipid.copy()
    
    ax.set_xlabel('Cluster Number', fontsize = 20, fontweight = 'bold')
    ax.set_ylabel('Number of Occurrences', fontsize = 20, fontweight = 'bold')
    
    ax.set_xticks(range(1,nCP.shape[0]+1))
    ax.tick_params(labelsize = 18)
    
    ax.legend(fontsize = 18)
    
    fig.tight_layout()

ION = ['No Ions', 'With Ions']
CP  = ['PC', 'PE', 'PG', 'PS', 'PC-PE', 'PC-PG', 'PC-PS', 'PE-PG', 'PE-PS', 'PG-PS']

with open('data_apl_ions.pickle', 'rb') as f:
    APL = pickle.load(f)
f.close()

with open('data_curvature_ions.pickle', 'rb') as f:
    CUR = pickle.load(f)
f.close()

with open('data_thickness_ions.pickle', 'rb') as f:
    THI = pickle.load(f)
f.close()

with open('data_tails_ions.pickle', 'rb') as f:
    TAI = pickle.load(f)
f.close()

labels = []

force = []
compo = []

lipids = {'PC':[], 'PE':[], 'PG':[], 'PS':[]}


Outliers = []

for idx, prop in enumerate(list([APL, CUR, THI, TAI])):
    
    kk = []
    
    for comp in range(prop.shape[1]):
        
        out    = np.zeros((prop.shape[0]), dtype = str)
        out[:] = str('N')
        
        data     = prop[:,comp]
        q25, q75 = np.percentile(data, [25, 75])
        iqr      = q75-q25
        
        ids_min = np.where(data < q25 - 1.5*iqr)[0]
        ids_max = np.where(data > q75 + 1.5*iqr)[0]
        
        ids = np.append(ids_min, ids_max)
        
        out[ids] = str('Y')
        
        kk.append(out)
    
    kk = np.array(kk)
    kk = kk.T
    
    Outliers.append(kk)

Outliers = np.array(Outliers)

dataframe = np.zeros((192*20,4))

outliers = {'APL':[], 'CUR':[], 'THI':[], 'TAI':[]}

n = 0
for apl, cur, thi, tai, out_apl, out_cur, out_thi, out_tai in zip(APL, CUR, THI, TAI, Outliers[0], Outliers[1], Outliers[2], Outliers[3]):
    
    l = apl.shape[0]
    
    dataframe[n*l:(n+1)*l,0] += apl
    dataframe[n*l:(n+1)*l,1] += cur
    dataframe[n*l:(n+1)*l,2] += thi
    dataframe[n*l:(n+1)*l,3] += tai
    
    for oapl, ocur, othi, otai in zip(out_apl, out_cur, out_thi, out_tai):
        outliers['APL'].append(oapl)
        outliers['CUR'].append(ocur)
        outliers['THI'].append(othi)
        outliers['TAI'].append(otai)
        
    for ff in ION:
        for cp in CP:
            
            labels.append(ff+' '+cp)
            force.append(ff)
            compo.append(cp)
            
            for lipid in lipids:
                
                if lipid in cp:
                    if lipid == cp:
                        lipids[lipid].append(1)
                    else:
                        lipids[lipid].append(0.5)
                else:
                    lipids[lipid].append(0)    
    
    n+=1


labels = np.array(labels)

df = pd.DataFrame(data = dataframe,
                  index = labels,
                  columns = ['APL', 'Curvature', 'Thickness', 'Tail Tilt'])

scaler = StandardScaler()
scaler.fit(df)
scaled = scaler.fit_transform(df)

df = pd.DataFrame(data = scaled,
                  index = labels,
                  columns =  ['APL', 'Curvature', 'Thickness', 'Tail Tilt'])

df['APL Outliers']       = outliers['APL']
df['Curvature Outliers'] = outliers['CUR']
df['Thickness Outliers'] = outliers['THI']
df['Tilt Outliers']      = outliers['TAI']

df['POPC']        = lipids['PC']
df['POPE']        = lipids['PE']
df['POPG']        = lipids['PG']
df['POPS']        = lipids['PS']

df['Ions']        = force
df['Composition'] = compo


id_AA = np.where(df['Ions'] == 'No Ions')[0]
id_CG = np.where(df['Ions'] == 'With Ions')[0]

dfAA = df.iloc[id_AA]
dfCG = df.iloc[id_CG]

clustered_total, nclust_total = Clusterize(df.copy(), 10)
clustered_AA, nclust_AA       = Clusterize(dfAA.copy(), 9)
clustered_CG, nclust_CG       = Clusterize(dfCG.copy(), 9)

Graphs(clustered_total, nclust_total, total = True, N = df.shape[0])

#Graphs(clustered_AA, nclust_AA, total = False, N = dfAA.shape[0])

#Graphs(clustered_CG, nclust_CG, total = False, N = dfCG.shape[0])
