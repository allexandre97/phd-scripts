import numpy as np
from MDAnalysis import *
import optparse
import fileinput
from math import *
import os

# the bilayer should not be broken; check the list of admited lipids, they should be in the table lipids of the supepmem DB

desc=''
phos = {'CG':'PO4', 'AA':'P'}

def getlipids(pdb):
    ress = []
    for res in pdb.residues:
        if res.resname != 'CHOL' and res.resname != 'NA' and res.resname != 'TIP3':
            if not res.resname in ress:
                ress.append(res.resname)

    return ress

def main() :
    parser = optparse.OptionParser(description=desc, version='%prog version 1.0')
    parser.add_option('-f', '--pdb_file', help='input pdb file', action='store')
    parser.add_option('-r', '--res', help='system resolution (AA or CG)', action='store')
    parser.set_defaults(pdb_file='equilibration3.gro')
    options, arguments = parser.parse_args()

#**************************************************************************************
    pdb_file=str(options.pdb_file)

    pdb=Universe(pdb_file)

    lpds = getlipids(pdb)

    p = ''
    for l in lpds:
        p += '(resname '+l+' and name '+phos[str(options.res)]+') '
        if l != lpds[-1]:
            p += 'or '

    P8=pdb.select_atoms(p) 

    zmax=np.max(P8.positions[:,2]) 
    zmin=np.min(P8.positions[:,2]) 
    zav=(zmax+zmin)*0.5

    lip = ''
    for l in lpds:
        lip += '(resname '+l+') '
        if l != lpds[-1]:
            lip += 'or '

    lipids=pdb.select_atoms(lip)

    resid=-1;resid=0; l_i=[]; leaflet=-1
    lowerndx=''; upperndx=''
    lowerpdb=''; uppperpdb=''
    for l in lipids:
        if l.resindex!=resid:
            resid=l.resindex
            if leaflet==0: 
                for i in l_i: lowerndx=lowerndx+str(i)+' '
                lowerndx=lowerndx+'\n'
            elif leaflet==1:
                for i in l_i: upperndx=upperndx+str(i)+' '
                upperndx=upperndx+'\n'

            l_i=[]
        if l.name==phos[options.res] and l.position[2]>=zmin and l.position[2]<zav:   leaflet=0
        elif l.name==phos[options.res] and l.position[2]<=zmax and l.position[2]>zav: leaflet=1
        l_i.append(l.id)

    if leaflet==0: 
        for i in l_i: lowerndx=lowerndx+str(i)+' '
    elif leaflet==1:
        for i in l_i: upperndx=upperndx+str(i)+' '


    f=open("leaflets.ndx", "w+")
    f.write("[ lower ]\n%s\n" % lowerndx)
    f.write("[ upper ]\n%s" % upperndx)
    f.close()

if __name__=="__main__" :
    main()

