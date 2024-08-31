# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 10:57:22 2023

@author: alexandre
"""

import MDAnalysis as mda

U = mda.Universe('TEMP_313/init.gro', 'whole.trr', in_memory = False)


for frame in U.trajectory[-100:]:
    
    print(U.atoms[0].velocity)
    print(U.atoms[0].force)


