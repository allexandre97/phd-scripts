#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 10:11:40 2024

@author: alexandre
"""

import numpy             as np
import pickle            as pk
import matplotlib.pyplot as plt

from scipy.stats         import mannwhitneyu
from matplotlib.gridspec import GridSpec as GS

sig_level = 0.001

MEMBRANES = ['POPC',
             'POPE',
             'POPG',
             'POPS']

PROPS = ['apl',
         'cur',
         'thi',
         'sn1a',
         'sn1b',
         'sn1c',
         'sn2a',
         'sn2b',
         'sn2c']


if __name__ == '__main__':
    
    DATA = {}
    
    for membrane in MEMBRANES:
        
        with open(f'./{membrane}/PRES.pk', 'rb') as f:
            PRES = pk.load(f)
        f.close()
    
        with open(f'./{membrane}/FREQ.pk', 'rb') as f:
            FREQ = pk.load(f)
        f.close()

        memdata = {}
        
        for prop in PROPS:
            
            prop_data = []
            
            for press in range(5, 55, 5):
                
                tmp = []
                
                for freq in range(5, 55, 5):
                    
                    master_path = f'{membrane}/{press}Bar_{freq}MHz/'
            
                    with open(f'{master_path}area_{prop}.pk', 'rb') as f:
                        data = pk.load(f)
                    f.close()
                    
                    tmp.append(data)
                
                prop_data.append(tmp)
            
            memdata[prop] = prop_data
        
        DATA[membrane] = memdata
        
    
    PVALS = np.zeros((4,4,10,10,9,3))
    
    for i, mem_i in enumerate(MEMBRANES):
        
        for j, mem_j in enumerate(MEMBRANES):
            
            for p in range(10):
                
                for f in range(10):
                    
                    for k, prop in enumerate(PROPS):
                        
                        print(i,j,p,f,k)
                        
                        data_i = DATA[mem_i][prop][p][f]
                        data_j = DATA[mem_j][prop][p][f]
                        
                        _, pval       = mannwhitneyu(data_i, data_j,
                                                     alternative = 'two-sided',
                                                     method = 'exact')
                        _, pval_less  = mannwhitneyu(data_i, data_j,
                                                     alternative = 'less',
                                                     method = 'exact')
                        _, pval_great = mannwhitneyu(data_i, data_j,
                                                     alternative = 'greater',
                                                     method = 'exact')
                        
                        PVALS[i,j,p,f,k] = np.array([pval, pval_less, pval_great])
                        
    
    with open('PVALS.pk', 'wb') as f:
        
        pk.dump(PVALS, f,
                protocol = 4)
    
    f.close() 
