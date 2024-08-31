# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 17:12:31 2023

@author: alexandre
"""

import numpy  as np
import pickle as pk
import MDAnalysis as mda
import MDAnalysis.transformations as trf

from MDAnalysis.analysis.leaflet import LeafletFinder as LF

intfmt = '{:>3}'

if __name__ == '__main__':
    print('alkfhalkfa')
    U = mda.Universe('./prod.tpr',
                     './prod.xtc', in_memory = False)
     
    DATA = []
    
    for f, frame in enumerate(U.trajectory):
        
        print(f'Done {intfmt.format(int(100*f/(U.trajectory.n_frames - 1)))} %',
              end = '\r')
        
        area = frame.dimensions[0]*frame.dimensions[1]

        DATA.append([frame.time, area])
    
    DATA = np.array(DATA)
    
    with open('area_total.pk', 'wb') as f:
        pk.dump(DATA, f,
                protocol = 4)
    f.close()
    
