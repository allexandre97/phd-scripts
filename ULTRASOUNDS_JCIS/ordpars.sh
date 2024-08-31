#!/bin/bash
#Set job requirements
#SBATCH --signal=TERM@5
#SBATCH -n 1
#SBATCH -c 50
#SBATCH -t 6:00:00
#SBATCH -p medium
#SBATCH --mem=246G

declare -a Lipids=( "POPC" "POPE" "POPG" "POPS" "POPC" "POPE" "POPC" "POPG" "POPC" "POPS" "POPE" "POPG" "POPE" "POPS" "POPG" "POPS" )
declare -a Folders=("PC" "PE" "PG" "PS" "PC-PE" "PC-PG" "PC-PS" "PE-PG" "PE-PS" "PG-PS")

Folder_0=`pwd`

n=0
for folder in ${Folders[@]}; do
	if [ $folder == "PC" ] || [ $folder == "PE" ] || [ $folder == "PG" ] || [ $folder == "PS" ]; then
		echo $folder
		echo ${Lipids[$n]}
		python ../../graph_ordp.py -l ${Lipids[$n]}
		n=$((n+1))
	else
		echo $folder
		echo ${Lipids[$n]}
		python ../../graph_ordp.py -l ${Lipids[$n]}
		n=$((n+1))
		echo $folder
		echo ${Lipids[$n]}
		python ../../graph_ordp.py -l ${Lipids[$n]}
		n=$((n+1))
	fi
	cd $Folder_0
done
