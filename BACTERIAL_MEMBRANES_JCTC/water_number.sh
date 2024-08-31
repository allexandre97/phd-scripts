#!/bin/bash


echo "Number of Water molecules inside the bilayer" > number_of_waters.txt
echo "The inner part of the bilayer is defined as the region between the glycerol lateral density peaks" >> number_of_waters.txt
echo "" >> number_of_waters.txt

for mem in ./POPC/ ./POPE/ ./POPG/ ./POPC-POPG_7-3/ ./POPE-POPG_1-3/ ./POPE-POPG_3-1/; do

  for ff in .; do

    cd ${mem}${ff}
    echo ${mem}${ff}

    : '
    cp ../../calc_thickness.py .
    python calc_thickness.py
    '

    cp ../water_dens.py .

    echo "$mem $ff" >> ../number_of_waters.txt
    python ./water_dens.py >> ../number_of_waters.txt
    echo "" >> ../number_of_waters.txt

    cd -
  done
done
