#!/bin/bash

#sed -i '/time_*/d' out*
#sed -i '/Application/d' out*.txt

for((i=0;i<100;i++))
do
	split -l 5000 out_indep$i.txt out_indep_split$i.txt
done

cat out_indep_split*.txtaa > time_iter_indep.txt
cat out_indep_split*.txtab > time_comp_indep.txt
cat out_indep_split*.txtac > time_comm_indep.txt
