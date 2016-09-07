#!/bin/bash
cd slabs/
for ((c=1;c<=$1;c++))
do
	sed -i '$ d' $c.grd
	sed -i '$ d' $c.grd
done
