#!/bin/bash

cd slabs

n=(*)
n=${#n[@]}

cd ../

x=''
i=0

file="grd.info"

if [ ! -f "$file" ] ; then
   touch $file
fi

echo $n > $file
