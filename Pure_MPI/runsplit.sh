#!/bin/bash

#sed -i '/time_*/d' out*
#sed -i '/Application/d' out*

for((i=0;i<100;i++))
do
	split -l 5000 out$i.txt out_split$i.txt
done

cat *.txtaa > time_iter.txt
cat *.txtab > time_comp.txt
cat *.txtac > time_comm.txt
