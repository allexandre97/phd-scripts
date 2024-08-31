# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 11:51:29 2023

@author: alexandre
"""

import MDAnalysis as mda
import numpy as np



if __name__ == '__main__':
    
    U = mda.Universe('100ns.gro', '2ns.xtc',
                     in_memory = False)

