#!/bin/bash

for rad in 040 050 060 080 100 120 140 160 180 200 220 250 #265;
do

	tail +13 '../physical_structure/'$rad'AU.txt' | while read z n_H Tgas Av Diff Tdust abun AVNH a;
	do
		echo "rad = $rad AU, elevation: $z AU"
		Temp=`python get_results_Dec2018.py -f 6.22e-5Msun_model.fits --coord $rad $z`
		echo "$Temp" | awk '{printf "%.3e ", $12}' | tr '\n' ' ' | cut -c11-170 >> T$rad.txt
		echo
	done
done
