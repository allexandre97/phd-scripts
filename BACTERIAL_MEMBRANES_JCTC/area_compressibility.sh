#!/bin/bash

module --force purge
module load cesga/2020 gcc/system openmpi/4.0.5_ft3 gromacs/2021.5

echo "Area Compressibility calculated with the area fluctuation method" > area_compressibility.txt
echo "" >> area_compressibility.txt

for mem in ./POPC/ ./POPE/ ./POPG/ ./POPC-POPG_7-3/ ./POPE-POPG_1-3/ ./POPE-POPG_3-1/; do

  for ff in .; do

    cd ${mem}${ff}
    echo ${mem}${ff}

    if [ ! -d EDR ]; then

	    mkdir EDR
	    mv *edr EDR/.

    fi

    rm EDR/total.edr
    gmx eneconv -f EDR/*edr -o EDR/total.edr -b 400000

    cp ../compressibility.py .

    echo "$mem $ff" >> ../area_compressibility.txt
    python ./compressibility.py -f EDR/total.edr >> ../area_compressibility.txt
    echo "" >> ../area_compressibility.txt

    cd -
  done
done
