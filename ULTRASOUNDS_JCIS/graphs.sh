#!/bin/bash

Perform_Analyses(){

        beg=${1:-0}
	end=${2:-500}
        processors=${3:-32}

	python ./Graphing/graph_apl.py $processors $beg $end
	python ./Graphing/graph_curv.py $processors $beg $end
	python ./Graphing/graph_flip.py $beg $end
	
	Lipids=`ls ordpa*pickle`
	for lipid in $lipids; do
		python ./Graphing/graph_ordp.py $lipid $processors $beg $end
	done
}


declare -a Lipids=( "POPC" "POPE" "POPG" "POPS" "POPC" "POPE" "POPC" "POPG" "POPC" "POPS" "POPE" "POPG" "POPE" "POPS" "POPG" "POPS" )
declare -a Folders=("PC" "PE" "PG" "PS" "PC-PE" "PC-PG" "PC-PS" "PE-PG" "PE-PS" "PG-PS")

Folder_0=`pwd`
n=0

TOTAL_JOBS=${#Folders[@]}
MAX_JOBS=${1:-1}
TOTAL_CORES=${2:-12}
BEG=${3:-0}
END=${4:-500}

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

module load cesga/2020 matplotlib/3.5.3-python-3.9.9

while [ $JOB_ID -le $TOTAL_JOBS ]; do

	if [ $n -eq $MAX_JOBS ]; then
		n=0
		wait
	fi

	folder=${Folders[$((JOB_ID - 1))]} &&
	
	cd $folder &&
	
	cp -r ../Graphing/ . &&

	echo `pwd`

	Perform_Analyses $BEG $END $CORES_PER_TASK &

	cd $Folder_0
	n=$((n + 1))
	JOB_ID=$((JOB_ID + 1))


done
