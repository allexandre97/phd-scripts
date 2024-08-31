#!/bin/bash
#Set job requirements
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -C clk
#SBATCH -t 24:00:00
#SBATCH -p short
#SBATCH --mem=32G

module --force purge
module load cesga/2020 gcccore/system python/2.7.18

declare -a MEMBRANES=(POPC POPE POPG POPC-POPG_7-3 POPE-POPG_1-3 POPE-POPG_3-1)
declare -a FORCEFIELDS=(.)
declare -a LIPIDS=(POPC POPE POPG POPC POPG POPE POPG POPE POPG)

for ff in ${FORCEFIELDS[@]}; do
    N=0
    for mem in ${MEMBRANES[@]}; do

        cd $mem/$ff

        cp ../orderpars.py .

        if [ "$mem" != "POPC" ] && [ "$mem" != "POPE" ] && [ "$mem" != "POPG" ]; then

	    echo "${LIPIDS[$N]}" "${LIPIDS[$(( N + 1 ))]}"
	    python2 orderpars.py -p prod.gro -x 400_500ns.xtc -i ../ITP/"${LIPIDS[$N]}"/Slipids/"${LIPIDS[$N]}".itp -a ignH -u "C29 C210" -o op_"${LIPIDS[$N]}".txt &&
	    python2 orderpars.py -p prod.gro -x 400_500ns.xtc -i ../ITP/"${LIPIDS[$(( N + 1 ))]}"/Slipids/"${LIPIDS[$(( N + 1 ))]}".itp -a ignH -u "C29 C210" -o op_"${LIPIDS[$(( N + 1 ))]}".txt &
            N=$(( N + 2 ))

        else
	    echo "${LIPIDS[$N]}"
	    python2 orderpars.py -p prod.gro -x 400_500ns.xtc -i ../ITP/"${LIPIDS[$N]}"/Slipids/"${LIPIDS[$N]}".itp -a ignH -u "C29 C210" -o op_"${LIPIDS[$N]}".txt &
            N=$(( N + 1 ))

        fi

        cd -
    done
    wait
done
