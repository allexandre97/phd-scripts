#!/bin/bash

module --force purge
module load cesga/2020 gcc/system openmpi/4.0.5_ft3 gromacs/2021.5

declare -a MEMBRANES=(POPC POPE POPG POPC-POPG_7-3 POPE-POPG_1-3 POPE-POPG_3-1)
declare -a FORCEFIELDS=(.)
declare -a LIPIDS=(POPC POPE POPG POPC POPG POPE POPG POPE POPG)


for ff in ${FORCEFIELDS[@]}; do
    N=0
    for mem in ${MEMBRANES[@]}; do

	cd $mem/$ff

	if [ "$mem" != "POPC" ] && [ "$mem" != "POPE" ] && [ "$mem" != "POPG" ]; then

	    { echo "\"P_&_${LIPIDS[$N]}\""; echo "\"P_&_${LIPIDS[$N]}\""; }  | gmx rdf -f prod.xtc -s prod.tpr -o rdf_${LIPIDS[$N]}_${LIPIDS[$N]}.xvg -n NDX/index_rdf.ndx -b 400000 -e 500000 -selrpos mol_com -seltype mol_com -rmax 2 &
	    { echo "\"P_&_${LIPIDS[$N]}\""; echo "\"P_&_${LIPIDS[$(( N + 1 ))]}\""; }  | gmx rdf -f prod.xtc -s prod.tpr -o rdf_${LIPIDS[$N]}_${LIPIDS[$(( N + 1 ))]}.xvg -n NDX/index_rdf.ndx -b 400000 -e 500000 -selrpos mol_com -seltype mol_com -rmax 2 &
	    { echo "\"P_&_${LIPIDS[$(( N + 1 ))]}\""; echo "\"P_&_${LIPIDS[$(( N + 1 ))]}\""; }  | gmx rdf -f prod.xtc -s prod.tpr -o rdf_${LIPIDS[$(( N + 1 ))]}_${LIPIDS[$(( N + 1 ))]}.xvg -n NDX/index_rdf.ndx -b 400000 -e 500000 -selrpos mol_com -seltype mol_com -rmax 2 &

	    N=$(( N + 2 ))

	else

	     { echo "\"P_&_${LIPIDS[$N]}\""; echo "\"P_&_${LIPIDS[$N]}\""; } | gmx rdf -f prod.xtc -s prod.tpr -o rdf_${LIPIDS[$N]}_${LIPIDS[$N]}.xvg -n NDX/index_rdf.ndx -b 400000 -e 500000 -selrpos mol_com -seltype mol_com -rmax 2 &
	    N=$(( N + 1 ))

	fi

	cd -

    done

done
