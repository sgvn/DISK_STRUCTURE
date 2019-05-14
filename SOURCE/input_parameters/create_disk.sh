#!/bin/bash

# has to be located in profile folder... in the usual folders.
n=2
nb_points=64
element1=H2
element2=CO
element3=N2
#element3=CS
#element4=CN
#element5=S+
#element6=H2CO
dir=gas
mkdir -p $dir
file1="$dir/r1.in"
file2="$dir/r2.in"
file3="$dir/r3.in"
autocm=1.49597e13

####CREATES STRUCTURE####
#######################################

for d in ./*AU/ ;
do 
	radau=`echo $d | sed "s/\.//g" | sed "s/\///g"`
	radius=`echo $d | sed 's/\///g' | sed 's/AU//g' | sed 's/\.//g' | awk '{printf "%5.3E\n", $1}'`
	#creates radius column
	for i in $(seq 1 $nb_points);
	do
		echo $radius >> $dir/radius.in
	done
	
	#creates z/H column
	for j in $(seq 0 $((nb_points-1)));
	do
		z=$(echo "4*(1-(2*$j/127))" | bc -l )
		printf "%.6E\n" $z >> $dir/zH.in
	done

	
	#creates physical structure columns
	tail -n +13 "$d"/1D_static.dat | awk '{print $1,$2,$3,$4}' >> $dir/physical.in


	#creates abundance column of element1
	tail -n 1 "$d"/ab/$element1.ab | tr -s ' ' '\n' | tail -n +3 >> $dir/ab1.in 
	
	#for elt in $(seq 1 5);
	#do
	#	tail -n 1 "$d"/ab/$element$elt.ab | tr -s ' ' '\n' | tail -n +3 >> plots/2Dplot_ab$elt.in 
	#	echo >>plots/2Dplot_ab$elt.in
	#done
	
	#create abundance column of element2
	tail -n 1 "$d"/ab/$element2.ab | tr -s ' ' '\n' | tail -n +3 >> $dir/ab2.in 
	
	#create abundance column of element3
	tail -n 1 "$d"/ab/$element3.ab | tr -s ' ' '\n' | tail -n +3 >> $dir/ab3.in 
	
	#create abundance column of element4
	#tail -n 1 "$d"/ab/$element4.ab | tr -s ' ' '\n' | tail -n +3 >> $dir/2Dplot_ab4.in 
	#echo >>plots/2Dplot_ab4.in
	
	#create abundance column of element5
	#tail -n 1 "$d"/ab/$element5.ab | tr -s ' ' '\n' | tail -n +3 >> $dir/2Dplot_ab5.in 
	#echo >>plots/2Dplot_ab5.in
	
	#create abundance column of element6
	#tail -n 1 "$d"/ab/$element6.ab | tr -s ' ' '\n' | tail -n +3 >> $dir/2Dplot_ab6.in 
	#echo >>$dir/2Dplot_ab6.in


	paste -d" " $dir/radius.in $dir/zH.in $dir/physical.in $dir/ab1.in $dir/ab2.in $dir/ab3.in >> $file1
	awk '{printf "%8.6e\n", $4*$7}' $file1 >> $dir/n1.in
	awk '{printf "%8.6e\n", $4*$8}' $file1 >> $dir/n2.in
	awk '{printf "%8.6e\n", $4*$9}' $file1 >> $dir/n3.in

	paste -d" " $file1 $dir/n1.in $dir/n2.in $dir/n3.in >> $file2

	z1=`awk '{print $3}' $file2 | tail -n +5 | head -n 1`
	z2=`awk '{print $3}' $file2 | tail -n +6 | head -n 1`
	dcm=`echo | awk 'END{printf "%8.6e", '$z1*$autocm-$z2*$autocm'}'`

	awk 'sum+=$10{printf "%8.6e", sum; printf "\n"}' $file2 >> $dir/1_$element1.in
	awk 'sum+=$11{printf "%8.6e", sum; printf "\n"}' $file2 >> $dir/1_$element2.in
	awk 'sum+=$12{printf "%8.6e", sum; printf "\n"}' $file2 >> $dir/1_$element3.in

	awk -F, '{$1=$1*a;print}' a=$dcm OFS=, $dir/1_$element1.in >> $dir/2_$element1.in
	awk -F, '{$1=$1*a;print}' a=$dcm OFS=, $dir/1_$element2.in >> $dir/2_$element2.in
	awk -F, '{$1=$1*a;print}' a=$dcm OFS=, $dir/1_$element3.in >> $dir/2_$element3.in

	awk '{printf "%8.6e\n", $1}' $dir/2_$element1.in >> $dir/c_col_$element1.in
	awk '{printf "%8.6e\n", $1}' $dir/2_$element2.in >> $dir/c_col_$element2.in
	awk '{printf "%8.6e\n", $1}' $dir/2_$element3.in >> $dir/c_col_$element3.in
	

	paste -d" " $file2 $dir/c_col_$element1.in $dir/c_col_$element2.in $dir/c_col_$element3.in >> $file3

	echo "#radius z/H z n_H T av ab($element1) ab($element2) ab($element3) n($element1) n($element2) n($element3) c_Col_$element1 c_Col_$element2 c_Col_$element3" >> $dir/$radau.txt       
	cat $file3 >> $dir/$radau.txt

	
	####CLEANING####
	################
	rm $file1 $file2 $file3 $dir/radius.in $dir/zH.in $dir/physical.in $dir/ab1.in $dir/ab2.in $dir/ab3.in $dir/n1.in $dir/n2.in $dir/n3.in 
	rm $dir/c_col_$element1.in $dir/c_col_$element2.in $dir/c_col_$element3.in
	rm $dir/1_$element1.in $dir/1_$element2.in $dir/1_$element3.in
	rm $dir/2_$element1.in $dir/2_$element2.in $dir/2_$element3.in
	################ 

done

echo "#radius z/H z n_H T av ab($element1) ab($element2) ab($element3) n($element1) n($element2) n($element3) c_Col_$element1 c_Col_$element2 c_Col_$element3" >> 2D_plot.txt
tail -n +2 $dir/*AU.txt | grep -v '=' >> 2D_plot.txt
 