#!/bin/bash

module --force purge
module load cesga/2020 gcc/system openmpi/4.0.5_ft3 gromacs/2021.5


calc_hbonds(){

  max_n=${#components[@]}
  max_n=$((max_n - 1))

  for i in `seq 0 $max_n`; do
    for j in `seq $i $max_n`; do
        c1=${components[$i]}
        c2=${components[$j]}
        if [ $c1 == "SOL" ] && [ $c2 == "SOL" ]; then
          continue
        fi
        { echo $c1; echo $c2; } | gmx hbond -f $mem/$ff/400_500ns.xtc -s $mem/$ff/prod.tpr \
                                            -n $mem/$ff/NDX/index_hbonds.ndx \
                                            -num $mem/$ff/HBND/hbond_${c1}_${c2}.xvg \
                                            -nthreads 2
    done
  done
}

declare -a MEMBRANES=(POPC POPE POPG POPC-POPG_7-3 POPE-POPG_1-3 POPE-POPG_3-1)
declare -a FORCEFIELDS=(.)
declare -a LIPIDS=(PC PE PG PC PG PE PG PE PG)
declare -a PARTS=(Headgroup Glycerol SOL)

for ff in ${FORCEFIELDS[@]}; do

  N=0

  for mem in ${MEMBRANES[@]}; do

    unset components

    if [ $mem == "POPC" ] || [ $mem == "POPE" ] || [ $mem == "POPG" ]; then

      declare -a components

      components[0]="Headgroup"_${LIPIDS[$N]}
      components[1]="Glycerol"_${LIPIDS[$N]}
      components[2]="SOL"
      N=$((N + 1))

      calc_hbonds &

    else

      declare -a components

      components[0]="Headgroup"_${LIPIDS[$N]}
      components[1]="Glycerol"_${LIPIDS[$N]}
      N=$((N+1))
      components[2]="Headgroup"_${LIPIDS[$N]}
      components[3]="Glycerol"_${LIPIDS[$N]}
      components[4]="SOL"
      N=$((N + 1))

      calc_hbonds &

    fi
  done
done

wait
