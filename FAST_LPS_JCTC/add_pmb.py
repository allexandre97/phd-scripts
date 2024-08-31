# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:41:25 2023

@author: alexandre
"""

import numpy      as np
import MDAnalysis as mda

from math                    import pi, sqrt, sin
from optparse                import OptionParser
from MDAnalysis              import transformations as trf
from scipy.spatial.transform import Rotation        as rot

parser = OptionParser(usage       = 'python wavesim.py -i <input build file> -o <output name>'
                                    'asdhkalkfjalkfj',
                      prog        = 'WaveSim',
                      description = 'This is a program used to simulate propagating waves through a medium.')

parser.add_option('-i', '--infile',
                  action = 'store', type = 'string',
                  default = 'membrane',
                  help = 'Input Filename')

parser.add_option('-o', '--out',
                  action = 'store', type = 'string',
                  default = 'membrane',
                  help = 'Output Filename')

parser.add_option('-n', '--nmol',
                  action = 'store', type = 'int',
                  default = 10,
                  help = 'Number of Molecules to insert')

options, args = parser.parse_args()

def PBCDist(P1, P2, L):
    d = P2 - P1
    d -= np.round(d/L)*L
    return np.linalg.norm(d)


def GenRandomPos(REGION, NMol, DIM, Cutoff):
    POS = []
    n = 0
    while n < NMol:
        accept = True
        px = np.random.uniform(REGION[0,0], REGION[0,1])
        py = np.random.uniform(REGION[1,0], REGION[1,1])
        pz = np.random.uniform(REGION[2,0], REGION[2,1])
        p = np.array([px, py, pz])
        for pos in POS:
            d = PBCDist(p, pos, DIM)
            if d < Cutoff:
                accept = False
            break
        if accept:
            POS.append(p)
            n += 1

    POS = np.array(POS)
    return POS

def MakeUniverse(N_Molecules, U0, P0,
                 resindex_offset = 0):

    POS0  = U0.atoms.positions.copy()
    POS0 -= U0.atoms.center_of_mass()

    n_atoms    = len(U0.atoms)    * N_Molecules
    n_residues = len(U0.residues) * N_Molecules

    resindices = []
    for n in range(P0.shape[0]):
        tmp = np.repeat(n, len(U0.atoms))
        for t in tmp:
            resindices.append(t)

    atomnames = [atom.name for atom in U0.atoms] * N_Molecules
    resnames  = [res.resname for res in U0.residues] * N_Molecules

    U = mda.Universe.empty(n_atoms,
                           n_residues=n_residues,
                           atom_resindex=resindices,
                           trajectory=True)

    U.add_TopologyAttr('resid', np.arange(resindex_offset, resindex_offset + n_residues))
    U.add_TopologyAttr('name', atomnames)
    U.add_TopologyAttr('resname', resnames)

    positions = []

    for p in P0:

        tmp = POS0.copy()

        anglex = np.random.uniform(0, 2*np.pi)
        angley = np.random.uniform(0, 2*np.pi)
        anglez = np.random.uniform(0, 2*np.pi)

        R = rot.from_euler('XYZ', [anglex, angley, anglez])

        tmp = R.apply(tmp)

        tmp[:, 0] += p[0]
        tmp[:, 1] += p[1]
        tmp[:, 2] += p[2]

        for t in tmp:
            positions.append(t)
    positions = np.array(positions)
    U.atoms.positions = positions

    return U

if __name__ == '__main__':

    U_MEM = mda.Universe(f'{options.infile}.pdb', in_memory = True)
    U_PMB = mda.Universe('PMB.pdb',               in_memory = True)

    DIMS = U_MEM.dimensions

    REMP = U_MEM.select_atoms('resname REMP')
    PPLP = U_MEM.select_atoms('resname POP* or resname CDL2')

    REMP_COM = REMP.center_of_mass()
    PPLP_COM = PPLP.center_of_mass()

    D = REMP_COM - PPLP_COM

    if D[2] > 0:

        ZMIN = REMP.positions[:,2].max()
        ZMAX = DIMS[2]

        REGION = np.array([[0, DIMS[0]],
                           [0, DIMS[1]],
                           [ZMIN*(1+0.1), ZMAX*(1-0.1)]])

    else:

        ZMIN = 0
        ZMAX = REMP.positions[:,2].min()

        REGION = np.array([[0, DIMS[0]],
                           [0, DIMS[1]],
                           [ZMIN*(1+0.1), ZMAX*(1-0.1)]])

    POS0 = GenRandomPos(REGION, options.nmol, DIMS[:3], 20)
    NEW_PMB = MakeUniverse(options.nmol, U_PMB, POS0)

    Merge = mda.Merge(U_MEM.atoms, NEW_PMB.atoms)

    transform = trf.boxdimensions.set_dimensions(DIMS)
    Merge.trajectory.add_transformations(transform)

    with mda.coordinates.PDB.PDBWriter(f'{options.out}.pdb') as writer:
        writer.write(Merge)
    writer.close()

    with open(f'{options.infile}_topol.top', 'a') as f:
        f.write(f'PMB\t{options.nmol}\n')
    f.close()
