#!/bin/bash
#Set job requirements
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -C clk
#SBATCH -t 08:00:00
#SBATCH -p short
#SBATCH --mem=16G

module --force purge
module load cesga/2020 gcc/system openmpi/4.0.5_ft3 gromacs/2021.5

declare -a MEMBRANES=(POPC POPE POPG POPC-POPG_7-3 POPE-POPG_1-3 POPE-POPG_3-1)
declare -a FORCEFIELDS=(.)
declare -a LIPIDS=(POPC POPE POG POPC POG POPE POG POPE POG)

for t0 in {400000..490000..10000}; do
    
    t1=$(( t0 + 10000 ))

    for ff in ${FORCEFIELDS[@]}; do
        N=0
        for mem in ${MEMBRANES[@]}; do

	    cd $mem/$ff

	    cp ../calc_densXY.py .

	    if [ "$mem" != "POPC" ] && [ "$mem" != "POPE" ] && [ "$mem" != "POPG" ]; then

	        python calc_densXY.py -x 400_500ns.xtc -s "${LIPIDS[$N]}_UP" -i $t0 -f $t1 -p . &&
	        python calc_densXY.py -x 400_500ns.xtc -s "${LIPIDS[$N]}_DN" -i $t0 -f $t1 -p . &&
	        
	        python calc_densXY.py -x 400_500ns.xtc -s "${LIPIDS[$(( N + 1 ))]}_UP" -i $t0 -f $t1 -p . &&
	        python calc_densXY.py -x 400_500ns.xtc -s "${LIPIDS[$(( N + 1 ))]}_DN" -i $t0 -f $t1 -p . &

	        N=$(( N + 2 ))

	    else

	        python calc_densXY.py -x 400_500ns.xtc -s "${LIPIDS[$N]}_UP" -i $t0 -f $t1 -p . &&
	        python calc_densXY.py -x 400_500ns.xtc -s "${LIPIDS[$N]}_DN" -i $t0 -f $t1 -p . &
	        
	        N=$(( N + 1 ))

	    fi

	    cd -
	done
	wait
    done
done
