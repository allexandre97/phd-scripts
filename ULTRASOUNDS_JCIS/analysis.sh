#!/bin/bash

Perform_Analyses(){

        initial_gro=${1:-initial.gro}
        traj_xtc=${2:-center.xtc}
        processors=${3:-32}
	membrane=$4

	echo "Performing Memsearch, $membrane Membrane"
	echo "Performing Memsearch, $membrane Membrane" &>analyses_log.txt
	echo &>>analyses_log.txt

        python ./Scripts/memsearch.py -p $initial_gro -x $traj_xtc -r CG -n $processors #&>>analyses_log.txt
	echo "Performing Surface Calculations, $membrane Membrane"
        echo "Performing Surface Calculations, $membrane Membrane" &>>analyses_log.txt
        echo &>>analyses_log.txt

        python ./Scripts/surface.py -k . -n $processors #&>>analyses_log.txt 
        
	echo "Performing Thickness, $membrane Membrane"
        echo "Performing Thickness, $membrane Membrane" &>>analyses_log.txt
        echo &>>analyses_log.txt

	python ./Scripts/thickness.py -k . -o thickness -n $processors #&>>analyses_log.txt
        
	echo "Performing Filp Flop, $membrane Membrane"
        echo "Performing Filp Flop, $membrane Membrane" &>>analyses_log.txt
        echo &>>analyses_log.txt

	python ./Scripts/flipflop.py -p $initial_gro -x $traj_xtc -r CG -k . -o flipflop -n $processors #&>>analyses_log.txt
        
	echo "Performing Area Per Lipid, $membrane Membrane"
        echo "Performing Area Per Lipid, $membrane Membrane" &>>analyses_log.txt
        echo &>>analyses_log.txt

	python ./Scripts/area_per_lipid.py -p $initial_gro -x $initial_xtc -r CG -o apl -n $processors #&>>analyses_log.txt
        
	echo "Performing Order Parameters, $membrane Membrane"
        echo "Performing Order Parameters, $membrane Membrane" &>>analyses_log.txt
        echo &>>analyses_log.txt

	python ./Scripts/orderpars_CoarseGrain.py -p $initial_gro -x $traj_xtc -k . -o ordpars -n $processors #&>>analyses_log.txt

}


declare -a Lipids=( "POPC" "POPE" "POPG" "POPS" "POPC" "POPE" "POPC" "POPG" "POPC" "POPS" "POPE" "POPG" "POPE" "POPS" "POPG" "POPS" )
declare -a Folders=("PC" "PE" "PG" "PS" "PC-PE" "PC-PG" "PC-PS" "PE-PG" "PE-PS" "PG-PS")

Folder_0=`pwd`
n=0

TOTAL_JOBS=${#Folders[@]}
MAX_JOBS=${1:-1}
TOTAL_CORES=${2:-12}

CORES_PER_TASK=`echo "import numpy as np; print(int(np.floor($TOTAL_CORES / $MAX_JOBS)))" | python`

echo
echo "#######################################################"
echo 
echo "STARTING ANALYSES:"
echo "------------------"
echo "--> Will Perform $TOTAL_JOBS Jobs"
echo
echo "--> Will Carry $MAX_JOBS Jobs in Parallel, using $CORES_PER_TASK Cores per Job"
echo 
echo "#######################################################"
echo

JOB_ID=1
n=0

while [ $JOB_ID -le $TOTAL_JOBS ]; do

	if [ $n -eq $MAX_JOBS ]; then
		n=0
		wait
	fi

	folder=${Folders[$((JOB_ID - 1))]} &&
	
	cd $folder &&
	
	cp -r ../Scripts/ . &&

	echo `pwd`

	Perform_Analyses initial.gro center.xtc $CORES_PER_TASK $folder &

	cd $Folder_0
	n=$((n + 1))
	JOB_ID=$((JOB_ID + 1))


done
