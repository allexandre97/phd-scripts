# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 13:53:55 2023

@author: alexandre
"""

import os
import numpy             as np
import pickle            as pk
import cblind            as cb
import matplotlib.pyplot as plt

from scipy.interpolate import CubicSpline
from scipy.integrate import simps

from statsmodels.tsa.stattools import acf

class Velocities:

    def __init__(self,
                 BOUND,
                 VELS):

        self.BOUND = BOUND
        self.VELS  = VELS

def Integrate(X, dX):

    out = [X[0]]

    for i, x in enumerate(X[1:]):

        out.append(out[i]+x)

    return np.array(out)*dX

with open('velocities.pk', 'rb') as f:
    VELS = pk.load(f)
f.close()

dt    = 20000/VELS.VELS.shape[0]

nlags = int(np.ceil(20/dt))

ACF = np.zeros(((VELS.VELS.shape[0]-1)*2+1))

for nmol in range(VELS.VELS.shape[2]):

    mol = VELS.VELS[:,0,nmol]

    for c in range(3):

        ACF += np.correlate(mol[:,c],
                            mol[:,c],
                            'full')


ACF /= np.sum(VELS.VELS[:,0,:,:]**2)

ACF = ACF[int(ACF.shape[0]/2):int(ACF.shape[0]/2)+nlags]

X = np.linspace(0, 20, nlags)

spline = CubicSpline(X, ACF)

Xnew = np.linspace(0, 20, 10*nlags)

print((1/VELS.VELS.shape[2])*simps(spline(Xnew), Xnew, dx = Xnew[1]-Xnew[0]))

fig, ax = plt.subplots()
ax2 = ax.twinx()
ax.scatter(X, ACF, c = 'royalblue')
ax.plot(Xnew, spline(Xnew), c = 'royalblue')
ax2.plot(Xnew, Integrate(spline(Xnew), Xnew[1]-Xnew[0]), c = 'firebrick')

