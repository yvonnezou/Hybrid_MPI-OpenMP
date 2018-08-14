#!/bin/bash

#sed -i '/time_*/d' out*
sed -i '/Application/d' out*.txt

for((i=0;i<100;i++))
do
	split -l 5000 outindep$i.txt out_splitindep$i.txt
done

cat out_splitindep*.txtaa > time_iter_indep.txt
cat out_splitindep*.txtab > time_comp_indep.txt
cat out_splitindep*.txtac > time_comm_indep.txt
