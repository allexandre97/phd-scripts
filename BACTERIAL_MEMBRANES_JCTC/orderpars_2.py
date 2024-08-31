#!/usr/bin/env python3

# Please, contact to Angel.Pineiro at usc.es for any question or comment
# Version May 9, 2019

import numpy as np
from MDAnalysis import Universe
import argparse
import fileinput
from math import *
import io

desc = 'Calculation of order parameters from a pdb file or a xtc trajectory. It works for AA and UA force fields. Explicit hydrogens can be ignored by using -a ignH. For UA force fields the unsaturated carbon atoms should be introduced using the -u flag. The atom list is taken from the topology itp file of the target molecule. Example: python2.7 opAAUA.py -p 1.pdb -i POPC.itp -r POPC -a ignH -u "C29 C210"'

# Examples:
# Single pdb file ignoring the explicit H atoms:    python2.7 opAAUA.py -p 1.pdb -i POPC.itp -r POPC -a ignH -u "C29 C210"
# Whole trajectory considering explicit atoms:      python2.7 opAAUA.py -p popc.pdb -i POPC.itp -r POPC -a AA -x popc.xtc
# Whole trajectory ignoring the explicit H atoms:   python2.7 opAAUA.py -p popc.pdb -i POPC.itp -r POPC -a ignH -u "C29 C210" -x popc.xtc

def permuta(string):
    aux = string[-1]
    for i in range(0, len(string) - 1, 1):
        aux = aux + string[i]
    return aux

def ReadITPfile(itp_file):
    control = -1
    atoms = []
    bonds = []
    for line in fileinput.input(itp_file):
        if line.strip():
            fields = line.split()
            if (fields[0] == '[') and (fields[2] == ']') and (fields[1] in ["atomtypes", "defaults"]):
                control = 0
            elif (fields[0] == '[') and (fields[2] == ']') and (fields[1] == "atoms"):
                control = 1
            elif (fields[0] == '[') and (fields[2] == ']') and (fields[1] == "bonds"):
                control = 2
            elif (fields[0] == '[') and (fields[2] == ']') and (fields[1] == "pairs"):
                control = 3
            elif (fields[0] == '[') and (fields[2] == ']') and (fields[1] == "angles"):
                control = 4
            # write to "atoms" and "bonds"
            if (control == 1) and (fields[0][0] != ';') and (fields[0] != '[') and (len(fields) > 1):
                atoms.append(fields)
            elif (control == 2) and (fields[0][0] != ';') and (fields[0] != '[') and (len(fields) > 1):
                bonds.append(fields)
    return [atoms, bonds]

def GetCatoms_array(itp, AA_UA, unsat):
    Catoms = GetCatoms(itp, AA_UA, unsat)
    C2atoms = GetC2atoms(Catoms)
    Carbon = ['C', 'c']
    aux = [catom for catom in Catoms if catom[0][0] in Carbon and len(catom) > 1]
    Catoms = aux
    aux = [c2atom for c2atom in C2atoms if c2atom[0][0] in Carbon]
    C2atoms = aux
    Catoms_per = Get_Catoms_per(Catoms)
    C2atoms_per = Get_Catoms_per(C2atoms)
    return (Catoms, Catoms_per, C2atoms, C2atoms_per)

def GetCatoms(itp, AA_UA, unsat):  # list of heavy atoms and the atoms to which they are bound, except Hydrogens if "-a ignH"
    atoms, bonds = itp
    Catoms = []
    Hydrogen = ['H', 'h']
    UA = ['ignH', 'IGNH']
    for a in range(len(atoms)):
        if atoms[a][1][0] not in Hydrogen:  # Assuming that type starting by H is a Hydrogen atom
            Catom_i = [atoms[a][4]]
            aid = atoms[a][0]
            for b in range(len(bonds)):
                aid2 = -1
                if bonds[b][0] == aid:
                    aid2 = bonds[b][1]
                elif bonds[b][1] == aid:
                    aid2 = bonds[b][0]
                if int(aid2) > 0:
                    for aa in range(len(atoms)):
                        if atoms[aa][0] == int(aid2):  # Assuming that type starting by H is a Hydrogen atom
                            if (AA_UA not in UA) and (atoms[aa][1][0] in Hydrogen):
                                Catom_i.append(atoms[aa][4])
                            elif (AA_UA in UA) and (atoms[aa][1][0] not in Hydrogen):
                                Catom_i.append(atoms[aa][4])
            Catoms.append(Catom_i)
    return Catoms

def Get_Catoms_per(Catoms):  # gromacs writes atomnames in the trajetory files with this permutation when the final character is a number
    numbers = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    Catoms_per = []
    for catom in Catoms:
        Catom_i_per = [permuta(name) if len(name) > 3 and name[-1] in numbers else name for name in catom]
        Catoms_per.append(Catom_i_per)
    return Catoms_per

def GetC2atoms(Catoms):  # list of second neighbours for C-atoms (only heavy atoms)
    C2atoms = []
    Hydrogen = ['H', 'h']
    for catom in Catoms:
        C2atom_i = [catom[0]]
        for catom2 in Catoms:
            if catom[0] in catom2[1:]:
                for i in range(1, len(catom2)):
                    if catom2[i] != catom[0] and catom2[i][0] not in Hydrogen:
                        C2atom_i.append(catom2[i])
        C2atoms.append(C2atom_i)
    return C2atoms

def GetOP(Catoms_array, trj, begin, end, skip, AA_UA, unsatlist, xtc_file, resname):
    op = []
    Catoms, Catoms_per, C2atoms, C2atoms_per = Catoms_array
    for c in range(len(Catoms)):
        if len(Catoms[c]) > 1:
            H = []
            X = []
            X2 = []
            op_c = []
            fr = -1
            Cname = Catoms[c][0]
            Cname_per = Catoms_per[c][0]
            Cnames = [Cname, Cname_per]
            print("\nAtom ", Cname, "\tRunning over frames ")
            if 1 == 2:#xtc_file == 'none':  # Single pdb file
                C = trj.select_atoms(f"resname {resname} and (name {Catoms[c][0]} or name {Catoms_per[c][0]})").positions
                for i in range(1, len(Catoms[c])):
                    if AA_UA == 'AA':
                        H.append(trj.select_atoms(f"resname {resname} and (name {Catoms[c][i]} or name {Catoms_per[c][i]})").positions)
                    else:
                        X.append(trj.select_atoms(f"resname {resname} and (name {Catoms[c][i]} or name {Catoms_per[c][i]})").positions)
                for i in range(1, len(C2atoms[c])):
                    X2.append(trj.select_atoms(f"resname {resname} and (name {C2atoms[c][i]} or name {C2atoms_per[c][i]})").positions)
                C = np.array(C)
                X = np.array(X)
                X2 = np.array(X2)
                H = np.array(H)
                if AA_UA != 'AA':
                    for j in range(len(X[0])):  # running over residues
                        op_c.append(GetOP_iUA(j, C, X, X2, Catoms, C2atoms, Cnames, unsatlist))
                else:
                    for j in range(len(H[0])):  # running over residues
                        op_c.append(GetOP_iAA(j, C, H[:, j, :], Catoms, Cnames))
                op.append(avop(Catoms[c], op_c))  # Average values
            else:  # Trajectory
                for frame in trj.trajectory[begin:end:skip]:  # Iterate over frames 
                    fr += 1
                    H = []
                    X = []
                    X2 = []
                    print(".", end="")
                    C = trj.select_atoms(f"resname {resname} and (name {Catoms[c][0]} or name {Catoms_per[c][0]})").positions
                    for i in range(1, len(Catoms[c])):
                        if AA_UA == 'AA':
                            H.append(trj.select_atoms(f"resname {resname} and (name {Catoms[c][i]} or name {Catoms_per[c][i]})").positions)
                        else:
                            X.append(trj.select_atoms(f"resname {resname} and (name {Catoms[c][i]} or name {Catoms_per[c][i]})").positions)
                    for i in range(1, len(C2atoms[c])):
                        X2.append(trj.select_atoms(f"resname {resname} and (name {C2atoms[c][i]} or name {C2atoms_per[c][i]})").positions)
                    if AA_UA != 'AA':
                        for j in range(len(X[0])):  # running over residues
                            op_c.append(GetOP_iUA(j, C, X, X2, Catoms, C2atoms, Cnames, unsatlist))
                    elif len(H) > 0:
                        H = np.array(H)
                        for j in range(len(H[0])):  # running over residues
                            op_c.append(GetOP_iAA(j, C, H[:, j, :], Catoms, Cnames))
                op.append(avop(Catoms[c], op_c))  # Average values
    return op

def GetOP_iUA(j, C, X, X2, Catoms, C2atoms, Cnames, unsatlist):
    C = np.array(C)
    X = np.array(X)
    X2 = np.array(X2)
    if Cnames[0] in unsatlist or Cnames[1] in unsatlist:
        H = GetH_unsat(Cnames, C[j], X[:, j, :], X2[:, j, :], Catoms, unsatlist)
    else:
        H = GetH(Cnames, C[j], X[:, j, :], X2[:, j, :], Catoms, C2atoms)
    op_i = [GetopH(C[j], h) for h in H]
    op_av = np.average(op_i)
    if len(H):
        op_i.append(op_av)
    if len(op_i) == 2:
        op_i.insert(1, 0)
    if len(op_i) == 3:
        op_i.insert(2, 0)
    if len(op_i) == 0:
        op_i = [0, 0, 0, 0]
    return op_i

def GetOP_iAA(j, C, H, Catoms, Cnames):
    C = np.array(C)
    H = np.array(H)
    op_i = [GetopH(C[j], h) for h in H]
    op_av = np.average(op_i)
    if len(H):
        op_i.append(op_av)
    if len(op_i) == 2:
        op_i.insert(1, 0)
    if len(op_i) == 3:
        op_i.insert(2, 0)
    if len(op_i) == 0:
        op_i = [0, 0, 0, 0]
    return op_i

def GetopH(C, H):
    d = np.linalg.norm(C - H)
    dz = H[2] - C[2]
    cos = dz / d
    return 0.5 * (3 * cos * cos - 1)

def GetH(Cnames, C, X, X2, Catoms, C2atoms):
    sin = 0.81664156
    cos = 0.57714519
    cos71 = np.cos(71 * np.pi / 180)
    sin71 = np.sin(71 * np.pi / 180)  # Constants
    V2 = C
    V1 = X[0]
    H = []
    if len(X) == 1:
        C2 = X2[0]  # Terminal C atom
        dv12 = V2 - V1
        dv12 = dv12 / np.linalg.norm(dv12)  # Unitary vector in the direction from X to C
        dv2C2 = V2 - C2
        dv2C2 = dv2C2 / np.linalg.norm(dv2C2)  # Unitary vector in the direction perpendicular to the C[i-2]-C[i-1]-C[i] plane
        Vp = np.cross(dv12, dv2C2)
        Vp = Vp / np.linalg.norm(Vp)  # Unitary vector perpendicular to dv12 (and to the C[i-2]-C[i-1]-C[i] plane)
        Vp2 = np.cross(dv12, Vp)
        Vp2 = Vp2 / np.linalg.norm(Vp2)  # Unitary vector perpendicular to dv12 and to Vp
        H3 = V2 + (cos71 * dv12 - sin71 * Vp2)  # H position in the plane formed by the C-chain and the vector joining C and X (109 deg with this vector)
        H.append(H3)
    elif len(X) == 2:
        V3 = X[1]
    if len(X) < 3:
        dv12 = V1 - V2
        dv23 = V2 - V3
        zV = V1 - V3  # Definition of vectors in the direction of X, Y, Z
        xV = np.cross(dv12, dv23)
        yV = np.cross(xV, zV)
        xV = xV / np.linalg.norm(xV)
        yV = yV / np.linalg.norm(yV)
        zV = zV / np.linalg.norm(zV)
        H1 = sin * xV + cos * yV
        H2 = -sin * xV + cos * yV
        H1 = H1 + V2
        H2 = H2 + V2  # Coordinates of the hydrogen atoms
        H.append(H2)
        H.append(H1)
    elif len(X) == 3:
        V1 = C
        V2 = X[0]
        V3 = X[1]
        V4 = X[2]
        dv12 = V2 - V1
        dv13 = V3 - V1
        dv14 = V4 - V1
        dv12 = dv12 / np.linalg.norm(dv12)
        dv13 = dv13 / np.linalg.norm(dv13)
        dv14 = dv14 / np.linalg.norm(dv14)
        H1 = -(dv12 + dv13 + dv14) + V1
        H.append(H1)
    return H

def GetH_unsat(Cnames, C, X, X2, Catoms, unsatlist):
    for a in range(len(Catoms)):
        if Catoms[a][0] in Cnames:
            C2name = Catoms[a][1]
            C3name = Catoms[a][2]
            C2name_per = permuta(C2name)
            C3name_per = permuta(C3name)
            C2 = X[0]
            C3 = X[1]
    if C2name in unsatlist or C2name_per in unsatlist:
        C2name, C3name = C3name, C2name
        C2, C3 = C3, C2
    sin = 0.5
    cos = float(sqrt(3.) / 2.)  # Constants
    V1 = C2
    V2 = C3
    V3 = C  # Coordinates of the three carbon atoms. The double bond is between V2 and V3
    dv12 = V1 - V2
    dv32 = zV = V3 - V2
    dv13 = V1 - V3  # Definition of vectors in the direction of X, Y, Z
    xV = np.cross(dv12, dv13)
    yV = np.cross(xV, zV)
    xV = xV / np.linalg.norm(xV)
    yV = yV / np.linalg.norm(yV)
    zV = zV / np.linalg.norm(zV)
    H = -cos * yV + sin * zV
    H = H + V3  # Coordinates of the hydrogen atom
    return [H]

def avop(Catom_i, op_c):
    op_c = np.array(op_c)
    c_i = [Catom_i[0]]
    for j in range(4):
        std = np.std(op_c[:, j])
        avg = np.average(op_c[:, j])
        c_i.append([avg, std, std / np.sqrt(len(op_c[:, j]))])
    return c_i

def write_results(op, outfile, Catoms):
    buf = io.StringIO()
    print("Atom_name  Hydrogen\tOP\t      STD\t   STDmean", file=buf)
    j = 0
    Carbon = ['C', 'c']
    for i in range(len(op)):
        if Catoms[i][0][0] in Carbon:
            j += 1
            if op[i][1][0]:
                print(f"{Catoms[i][0]:7s}\t{'HR':8s}  {op[i][1][0]:10.5f}\t{op[i][1][1]:10.5f}\t{op[i][1][2]:10.5f}", file=buf)
                if op[i][2][0]:
                    print(f"{'':7s}\t{'HS':8s}  {op[i][2][0]:10.5f}\t{op[i][2][1]:10.5f}\t{op[i][2][2]:10.5f}", file=buf)
                if op[i][3][0]:
                    print(f"{'':7s}\t{'HT':8s}  {op[i][3][0]:10.5f}\t{op[i][3][1]:10.5f}\t{op[i][3][2]:10.5f}", file=buf)
                print(f"{'':7s}\t{'AVG':8s}  {op[i][4][0]:10.5f}\t{op[i][4][1]:10.5f}\t{op[i][4][2]:10.5f}\n", file=buf)
    print(buf.getvalue())
    with open(outfile, "w") as fout:
        fout.write(buf.getvalue())

def main():
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-p', '--pdb_file', help='input pdb file', required=True)
    parser.add_argument('-x', '--xtc_file', help='input xtc file', default='none')
    parser.add_argument('-i', '--itp_file', help='input itp file', required=True)
    parser.add_argument('-b', '--begin', help='begin: first frame to be analyzed', type=int, default=0)
    parser.add_argument('-e', '--end', help='end: last frame to be analyzed', type=int, default=-1)
    parser.add_argument('-s', '--skip', help='skip: only write every nr-th frame', type=int, default=1)
    parser.add_argument('-a', '--AA_ignH', help='AA for All-Atoms and ignH to ignore all explicit H atoms', default='ignH')
    parser.add_argument('-u', '--unsat', help='list of unsaturated carbon atoms, only if "-a ignH"', default='')
    parser.add_argument('-o', '--outfile', help='output file name', default='out.op')
    args = parser.parse_args()

    pdb_file = str(args.pdb_file)
    xtc_file = str(args.xtc_file)
    itp_file = str(args.itp_file)
    begin = int(args.begin)
    end = int(args.end)
    skip = int(args.skip)
    outfile = str(args.outfile)
    AA_UA = str(args.AA_ignH)
    unsat = str(args.unsat)
    unsatlist = unsat.split()


    if xtc_file == 'none':
        trj = Universe(pdb_file)
    else:
        trj = Universe(pdb_file, xtc_file)
    if end == -1:
        end = trj.trajectory.n_frames

    itp = ReadITPfile(itp_file)
    resname = itp[0][0][3]
    Catoms_array = GetCatoms_array(itp, AA_UA, unsat)
    op = GetOP(Catoms_array, trj, begin, end, skip, AA_UA, unsatlist, xtc_file, resname)

    #print(f"\n\n\nAveraging order parameters for a total of {trj.residues.n_residues} lipids and {trj.trajectory.n_frames} frames\n")
    write_results(op, outfile, Catoms_array[0])

if __name__ == "__main__":
    main()
