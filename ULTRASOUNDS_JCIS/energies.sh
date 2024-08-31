#!/bin/bash

module load cesga/2020 gcc/system openmpi/4.0.5_ft3 gromacs/2021.5

for folder in PC PE PG PS PC-PE PC-PG PC-PS PE-PG PE-PS PG-PS; do
	cd $folder
	echo "Box-X
	Box-Y
	Box-Z
	Vir-XX
	Vir-XY
	Vir-XZ
	Vir-YX
	Vir-YY
	Vir-YZ
	Vir-ZX
	Vir-ZY
	Vir-ZZ
	Pres-XX
	Pres-XY
	Pres-XZ
	Pres-YX
	Pres-YY
	Pres-YZ
	Pres-ZX
	Pres-ZY
	Pres-ZZ
	T-MEMBRANE" | gmx energy -f complete.edr -o properties.xvg
	cd ..
done

