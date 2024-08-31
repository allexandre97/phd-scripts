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
    
    U = mda.Universe('./prod.gro',
                     './prod.xtc', in_memory = False)
    #U.trajectory.add_transformations(trf.unwrap(U.atoms))
    PO4   = U.select_atoms('name P*', periodic = True, updating = True)
    Leafs = LF(U, U.select_atoms('name P*', updating  = True), pbc = True) 

    Leaflet_A = Leafs.groups(0)
    Leaflet_B = Leafs.groups(1)
    
    DATA = []
    
    for f, frame in enumerate(U.trajectory):
        
        print(f'Done {intfmt.format(int(100*f/(U.trajectory.n_frames - 1)))} %',
              end = '\r')

        d = Leaflet_A.center_of_mass() - Leaflet_B.center_of_mass()
        
        DATA.append(np.linalg.norm(d))
    
    DATA = np.array(DATA)
    
    with open('thickness_total.pk', 'wb') as f:
        pk.dump(DATA, f,
                protocol = 4)
    f.close()
   
