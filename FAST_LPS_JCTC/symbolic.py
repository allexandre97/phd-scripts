# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 16:33:05 2023

@author: alexandre
"""

import sympy as sp
from sympy import Matrix, print_latex

a1b1, a1b2, a2b1, a2b2, a3b1, a3b2 = sp.symbols('p(a1b1) p(a1b2) p(a2b1) p(a2b2) p(a3b1) p(a3b2)')


M = sp.Matrix([[a1b1, a1b2],
               [a2b1, a2b2],
               [a3b1, a3b2]])

MT = M.T

M_dot_MT = sp.Matrix([[M.row(0).dot(MT.col(0)), M.row(0).dot(MT.col(1)), M.row(0).dot(MT.col(2)) ],
                      [M.row(1).dot(MT.col(0)), M.row(1).dot(MT.col(1)), M.row(1).dot(MT.col(2)) ],
                      [M.row(2).dot(MT.col(0)), M.row(2).dot(MT.col(1)), M.row(2).dot(MT.col(2)) ]])

print_latex(M_dot_MT)
