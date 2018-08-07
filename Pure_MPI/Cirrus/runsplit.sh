#!/bin/bash

#sed -i '/time_*/d' out*
#sed -i '/Application/d' out*

#for((i=0;i<100;i++))
#do
#	split -l 5000 out_indep$i.txt out_indep_split$i.txt
#done

cat out_split*.txtaa > time_iter.txt
cat out_split*.txtab > time_comp.txt
cat out_split*.txtac > time_comm.txt
