#!/bin/bash
cd tavolette/
for ((c=1;c<=$1;c++))
do
	x=''
	i=0
	while read line 
	do            
		IFS=' ' read -a array <<< "$line" 
		x[i]=${array[0]}
		((i++))
	done < $c.grd
	mv $c.grd ${x[2]}_${x[3]}.grd
done

