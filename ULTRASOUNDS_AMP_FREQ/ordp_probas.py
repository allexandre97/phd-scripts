#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 12:32:15 2024

@author: alexandre
"""

import numpy             as np
import pickle            as pk
import matplotlib.pyplot as plt

with open('POPC/35Bar_15MHz/ordp_POPC_mean.pickle', 'rb') as f:
    
    DATA = pk.load(f)

f.close()

with open('POPC/ShortEq/ordp_POPC_mean.pickle', 'rb') as f:
    
    DATA0 = pk.load(f)

f.close()


DATA0_mean = np.mean(DATA0, axis = (0, 2))[None,:,None,:]

DATA_norm  = DATA / DATA0_mean - 1

DATA_mean  = np.mean(DATA_norm, axis = 2)

#DATA  = DATA/DATA0 - 1

fig, ax = plt.subplots(figsize = (7,7))

ax.plot(DATA_mean[:,0,1])

ax.plot(DATA_mean[:,0,2])

