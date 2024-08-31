#!/usr/bin/bash

DensGroup=$1

module --force purge
module load cesga/2020 mdanalysis/2.1.0

for t0 in {400000..495000..5000}; do
  t1=$((t0 + 5000))
  for mem in POPC POPE POPG POPC-POPG_7-3 POPE-POPG_1-3 POPE-POPG_3-1; do
    for ff in .; do
      python ./calc_density.py -n $mem/$ff/NDX/index_dns.ndx -t $mem/$ff/prod.tpr \
                               -x $mem/$ff/400_500ns.xtc -s "Headgroup" -i $t0 -f $t1 \
                               -p $mem/$ff/ -o DensBootstrap/dens_hgp_${t0}_${t1}.pk &
      python ./calc_density.py -n $mem/$ff/NDX/index_dns.ndx -t $mem/$ff/prod.tpr \
                               -x $mem/$ff/400_500ns.xtc -s "Glycerol" -i $t0 -f $t1 \
                               -p $mem/$ff/ -o DensBootstrap/dens_glc_${t0}_${t1}.pk &
      python ./calc_density.py -n $mem/$ff/NDX/index_dns.ndx -t $mem/$ff/prod.tpr \
                               -x $mem/$ff/400_500ns.xtc -s "PO4" -i $t0 -f $t1 \
                               -p $mem/$ff/ -o DensBootstrap/dens_po4_${t0}_${t1}.pk &
      python ./calc_density.py -n $mem/$ff/NDX/index_dns.ndx -t $mem/$ff/prod.tpr \
                               -x $mem/$ff/400_500ns.xtc -s "Sn1" -i $t0 -f $t1 \
                               -p $mem/$ff/ -o DensBootstrap/dens_sn1_${t0}_${t1}.pk &
      python ./calc_density.py -n $mem/$ff/NDX/index_dns.ndx -t $mem/$ff/prod.tpr \
                               -x $mem/$ff/400_500ns.xtc -s "Sn2" -i $t0 -f $t1 \
                               -p $mem/$ff/ -o DensBootstrap/dens_sn2_${t0}_${t1}.pk &
      python ./calc_density.py -n $mem/$ff/NDX/index_dns.ndx -t $mem/$ff/prod.tpr \
                               -x $mem/$ff/400_500ns.xtc -s "SOL" -i $t0 -f $t1 \
                               -p $mem/$ff/ -o DensBootstrap/dens_sol_${t0}_${t1}.pk &
    done
  done
  wait
done
