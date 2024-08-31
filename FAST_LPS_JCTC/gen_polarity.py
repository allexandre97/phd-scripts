# mypy: ignore-errors

import numpy      as np
import MDAnalysis as mda

def GetPolarity(itpfile):

    with open(itpfile, 'r') as f:

        lines = f.readlines()

        POLARITY_ARRAY = []

        for l1, line in enumerate(lines):
            words = line.split()
            if len(words) > 1:
                if words[1] == 'atoms':
                    for l2, line2 in enumerate(lines[l1+1:]):
                        atomfields = line2.split()
                        if len(atomfields) == 0:
                            break
                        if atomfields[0] == ';' or atomfields[0][0] == ';':
                            continue
                        bead = atomfields[1]
                        if bead[0]=='C' or bead[0] == 'P':
                            pol = POLARITY[bead]
                        elif bead[0]=='S':
                            if bead[1] != 'N':
                                pol = POLARITY[bead[1:]]
                            else:
                                pol = POLARITY['N']
                        elif bead[0] == 'N':
                            pol = POLARITY['N']
                        elif bead[0] == 'Q':
                            pol = POLARITY['Q']

                        POLARITY_ARRAY.append(pol)

    return np.array(POLARITY_ARRAY)

PMB  = mda.Universe('./PMB.pdb', in_memory = True)
REMP = mda.Universe('./REMP.pdb', in_memory = True)

POLARITY = {'Q':0,  'N':0.5,
            'C5':1, 'C4':0.95, 'C3':0.90, 'C2':0.85, 'C1':0.80,
            'P5':0, 'P4':0.05, 'P3':0.10, 'P2':0.15, 'P1':0.20}

POL_PMB  = GetPolarity('./toppar/PMB.itp')
POL_REMP = GetPolarity('./toppar/REMP.itp')

PMB.add_TopologyAttr('tempfactors', POL_PMB)
REMP.add_TopologyAttr('tempfactors', POL_REMP)

with mda.coordinates.PDB.PDBWriter(f'POL_PMB.pdb') as writer:
    writer.write(PMB)
writer.close()

with mda.coordinates.PDB.PDBWriter(f'POL_REMP.pdb') as writer:
    writer.write(REMP)
writer.close()
