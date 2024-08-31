#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 12:12:25 2024

@author: alexandre
"""

import pickle as pk


with open('POPC/DATA.pk', 'rb') as f:
    
    data = pk.load(f)
    
f.close()

print(data[...,-1])
