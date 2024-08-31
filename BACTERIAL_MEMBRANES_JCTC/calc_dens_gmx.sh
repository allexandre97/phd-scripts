#!/bin/bash

module load cesga/2020 gcc/system openmpi/4.0.5_ft3 gromacs/2021.5

declare -a MEMBRANES=(POPC POPE POPG POPC-POPG_7-3 POPE-POPG_1-3 POPE-POPG_3-1)
declare -a FORCEFIELDS=(.)
declare -a LIPIDS=(PC PE PG PC PG PE PG PE PG)

for ff in ${FORCEFIELDS[@]}; do

  N=0

  for mem in ${MEMBRANES[@]}; do

    cd $mem/$ff

    if [ $mem == "POPC" ] || [ $mem == "POPE" ] || [ $mem == "POPG" ]; then

     { echo "PO4"; echo "Headgroup_${LIPIDS[$N]}"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_hgp_${LIPIDS[$N]}.xvg &
     { echo "PO4"; echo "Glycerol_${LIPIDS[$N]}"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_glc_${LIPIDS[$N]}.xvg &
     { echo "PO4"; echo "Sn1_${LIPIDS[$N]}"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_sn1_${LIPIDS[$N]}.xvg &
     { echo "PO4"; echo "Sn2_${LIPIDS[$N]}"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_sn2_${LIPIDS[$N]}.xvg &
     { echo "PO4"; echo "SOL"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_sol.xvg &
      N=$(( N+1 ))

    else

    { echo "PO4"; echo "Headgroup_${LIPIDS[$N]}"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_hgp_${LIPIDS[$N]}.xvg &
    { echo "PO4"; echo "Glycerol_${LIPIDS[$N]}"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_glc_${LIPIDS[$N]}.xvg &
    { echo "PO4"; echo "Sn1_${LIPIDS[$N]}"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_sn1_${LIPIDS[$N]}.xvg &
    { echo "PO4"; echo "Sn2_${LIPIDS[$N]}"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_sn2_${LIPIDS[$N]}.xvg &

      N=$(( N+1 ))

     { echo "PO4"; echo "Headgroup_${LIPIDS[$N]}"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_hgp_${LIPIDS[$N]}.xvg &
     { echo "PO4"; echo "Glycerol_${LIPIDS[$N]}"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_glc_${LIPIDS[$N]}.xvg &
     { echo "PO4"; echo "Sn1_${LIPIDS[$N]}"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_sn1_${LIPIDS[$N]}.xvg &
     { echo "PO4"; echo "Sn2_${LIPIDS[$N]}"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_sn2_${LIPIDS[$N]}.xvg &

     { echo "PO4"; echo "SOL"; } | gmx density -f 400_500ns.xtc -n NDX/index_dns.ndx -s prod.tpr -center -o DNS/latdens_sol.xvg &
      N=$(( N+1 ))

    fi

    cd -

  done

done

wait
