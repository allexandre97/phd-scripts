#!/bin/bash

dt=$1

declare -a MEMBRANES=(POPC POPE POPG POPC-POPG_7-3 POPE-POPG_1-3 POPE-POPG_3-1)
declare -a FORCEFIELDS=(.)
declare -a LIPIDS=(POPC POPE POG POPC POG POPE POG POPE POG)


for ff in ${FORCEFIELDS[@]}; do

    N=0

    if [ $ff == "Charmm" ] || [ $ff == "Slipids" ] || [ $ff == "." ]; then
	atom_name="P"
    else
	atom_name="P8"
    fi


    for mem in ${MEMBRANES[@]}; do

	cd $mem/$ff
	cp ../diff.py .

	
	if [ $mem == "POPC" ] || [ $mem == "POPE" ] || [ $mem == "POPG" ]; then

	    if [[ ( $ff == "GromosCKP"  ||  $ff == "GromosH2Q" )  &&  ${LIPIDS[$N]} == "POPG" ]]; then
	        residue="LPOG"
	    elif [ $ff == "Slipids" ] && [ ${LIPIDS[$N]} == "POPG" ]; then
	        residue="POG"
	    else
	        residue=${LIPIDS[$N]}
	    fi


	    python diff.py -p prod.gro -x 400_500ns.xtc -o dc_${LIPIDS[$N]}_all_${dt}ns.txt -d dc_${LIPIDS[$N]}_ind_${dt}ns.txt -t $dt -r ${residue} -a $atom_name &
	    N=$((N + 1))

	else

	    if [[ ( $ff == "GromosCKP"  ||  $ff == "GromosH2Q" )  &&  ${LIPIDS[$N]} == "POPG" ]]; then
	        residue="LPOG"
	    elif [ $ff == "Slipids" ] && [ ${LIPIDS[$N]} == "POPG" ]; then
	        residue="POG"
	    else
	        residue=${LIPIDS[$N]}
	    fi


	    python diff.py -p prod.gro -x 400_500ns.xtc -o dc_${LIPIDS[$N]}_all_${dt}ns.txt -d dc_${LIPIDS[$N]}_ind_${dt}ns.txt -t $dt -r ${residue} -a $atom_name &
            N=$((N + 1))

	    if [[ ( $ff == "GromosCKP"  ||  $ff == "GromosH2Q" )  &&  ${LIPIDS[$N]} == "POPG" ]]; then
	        residue="LPOG"
	    elif [ $ff == "Slipids" ] && [ ${LIPIDS[$N]} == "POPG" ]; then
	        residue="POG"
	    else
	        residue=${LIPIDS[$N]}
	    fi

	    python diff.py -p prod.gro -x 400_500ns.xtc -o dc_${LIPIDS[$N]}_all_${dt}ns.txt -d dc_${LIPIDS[$N]}_ind_${dt}ns.txt -t $dt -r ${residue} -a $atom_name &
            N=$((N + 1))
	fi

	cd -
    
    done
done 
