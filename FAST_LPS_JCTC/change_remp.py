# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 09:31:34 2023

@author: alexandre
"""

import numpy             as np
import MDAnalysis        as mda
import matplotlib.pyplot as plt

from MDAnalysis import transformations as trf

LPS_N_ATOMS = 42

def gen_COM_positions(AtomGroup, Atoms_Per_Molecule):

    PCOM = []

    N_Molecules = int(len(AtomGroup)/Atoms_Per_Molecule)

    for n in range(N_Molecules):
        tmp = AtomGroup[n*Atoms_Per_Molecule:(n+1)*Atoms_Per_Molecule]
        PCOM.append(tmp.center_of_mass())
    PCOM = np.array(PCOM)
    return PCOM


def MakeEqual(AtomGroup1, AtomGroup2, diff):

    allatoms = [atom for atom in AtomGroup1]
    for atom in AtomGroup2:
        allatoms.append(atom)

    #diff = int((len(AtomGroup1) - len(AtomGroup2))/2)

    new2 = [atom for atom in AtomGroup2]
    rand_atoms = np.random.choice(AtomGroup1, diff, replace = False)
    for atom in rand_atoms:
        new2.append(atom)
    new1 = [atom for atom in allatoms if not atom in new2]

    AtomGroup1_new = mda.AtomGroup(new1)
    AtomGroup2_new = mda.AtomGroup(new2)

    return AtomGroup1_new, AtomGroup2_new


def MakeUniverse(AtomGroup,
                 RESNAME,
                 NAME = None):

    n_atoms    = len(AtomGroup)
    n_residues = len(np.unique(AtomGroup.resids))

    resindices = np.repeat(np.arange(n_residues), int(n_atoms/n_residues))

    NewUniverse = mda.Universe.empty(n_atoms,
                                     n_residues    = n_residues,
                                     atom_resindex = resindices,
                                     trajectory    = True)

    resnames = [RESNAME for res in AtomGroup.resnames][:n_residues]

    if NAME == None:
        NewUniverse.add_TopologyAttr('name', AtomGroup.names)
    else:
        NewUniverse.add_TopologyAttr('name', [NAME]*n_atoms)

    NewUniverse.add_TopologyAttr('resname', resnames)
    NewUniverse.add_TopologyAttr('resid', np.arange(n_residues))
    NewUniverse.atoms.positions = AtomGroup.positions

    return NewUniverse


if __name__ == '__main__':

    U = mda.Universe('original.gro',
                     in_memory = True)

    REMP = U.select_atoms('resname REMP', periodic = True)


    REMP_COM = gen_COM_positions(REMP, LPS_N_ATOMS)
    ids_all  = np.arange(len(REMP_COM))


    ids_nonchanged = np.random.choice(ids_all, int(len(ids_all)/2), replace = False)
    ids_changed    = np.array([i for i in ids_all
                               if not i in ids_nonchanged])

    REMP_CHANGE    = mda.AtomGroup([atom for atom in REMP
                                    if atom.resid-1 in ids_changed])
    REMP_NONCHANGE = mda.AtomGroup([atom for atom in REMP
                                    if atom.resid-1 in ids_nonchanged])


    CALCIUM_CHANGE    = U.select_atoms('(resname ION and name CA) and around 4.9 group remp_change',
                                       remp_change = REMP_CHANGE, periodic = True)
    CALCIUM_NONCHANGE = U.select_atoms('(resname ION and name CA) and not around 4.9 group remp_change',
                                       remp_change = REMP_CHANGE, periodic = True)


    if len(CALCIUM_CHANGE) < len(CALCIUM_NONCHANGE):

        diff = int((len(CALCIUM_NONCHANGE)- len(CALCIUM_CHANGE))/2)

        CALCIUM_NONCHANGE, CALCIUM_CHANGE = MakeEqual(CALCIUM_NONCHANGE,
                                                      CALCIUM_CHANGE,
                                                      diff)

    elif len(CALCIUM_CHANGE) > len(CALCIUM_CHANGE):

        diff = int((len(CALCIUM_CHANGE)- len(CALCIUM_CHANGE))/2)

        CALCIUM_CHANGE, CALCIUM_NONCHANGE = MakeEqual(CALCIUM_CHANGE,
                                                      CALCIUM_NONCHANGE,
                                                      diff)

    assert(len(CALCIUM_CHANGE) == len(CALCIUM_NONCHANGE))


    #CALCIUM_NONCHANGE = U.select_atoms('resname ION and name CA')


    REMP_CHANGE   = MakeUniverse(REMP_CHANGE, 'FEMP')


    SODIUM_CHANGE = MakeUniverse(CALCIUM_CHANGE, 'ION', 'NA')

    POPE = U.select_atoms('resname POPE')
    POPG = U.select_atoms('resname POPG')
    CDL2 = U.select_atoms('resname CDL2')
    PMB1 = U.select_atoms('resname PMB1')
    WATR = U.select_atoms('resname W')
    SODM = U.select_atoms('resname ION and name NA')
    CHLR = U.select_atoms('resname ION and name CL')

    Merge = mda.Merge(REMP_NONCHANGE, REMP_CHANGE.atoms)
    Merge = mda.Merge(Merge.atoms, POPE)
    Merge = mda.Merge(Merge.atoms, POPG)
    Merge = mda.Merge(Merge.atoms, CDL2)
    Merge = mda.Merge(Merge.atoms, PMB1)
    Merge = mda.Merge(Merge.atoms, WATR)
    Merge = mda.Merge(Merge.atoms, CALCIUM_NONCHANGE)
    Merge = mda.Merge(Merge.atoms, SODIUM_CHANGE.atoms)
    Merge = mda.Merge(Merge.atoms, SODM)
    Merge = mda.Merge(Merge.atoms, CHLR)

    dim = U.dimensions
    transform = trf.boxdimensions.set_dimensions(dim)
    Merge.trajectory.add_transformations(transform)

    with mda.coordinates.PDB.PDBWriter('HalfFast/membrane.pdb') as writer:
        writer.write(Merge)
    writer.close()

    with open('topol_0.top', 'r') as f0:
        lines = f0.readlines()
        with open('HalfFast/system.top', 'w') as f1:
            for line in lines:
                f1.write(line)
            f1.write('\n')

            N_CA = len(Merge.select_atoms('resname ION and name CA'))
            N_NA = len(Merge.select_atoms('resname ION and name NA'))
            N_CL = len(Merge.select_atoms('resname ION and name CL'))

            f1.write(f'REMP\t{len(REMP_NONCHANGE.residues)}\n')
            f1.write(f'FEMP\t{len(REMP_CHANGE.residues)}\n')
            f1.write(f'POPE\t{len(POPE.residues)}\n')
            f1.write(f'POPG\t{len(POPG.residues)}\n')
            f1.write(f'CDL2\t{len(CDL2.residues)}\n')
            f1.write(f'PMB\t{len(PMB1.residues)}\n')
            f1.write(f'W\t{len(WATR.residues)}\n')
            f1.write(f'CA\t{N_CA}\n')
            f1.write(f'NA\t{N_NA}\n')
            f1.write(f'CL\t{N_CL}')
        f1.close()
    f0.close()


