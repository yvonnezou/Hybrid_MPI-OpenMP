#This file is for data processing including cleaning, splitting and gathering the data. 
#The data after processing will be stored in time_iter.txt, time_comp.txt and time_comm.txt.
#They store the runtime of each iteration, the computation time in each iteration and the communication time in each iteration respectively.

#!/bin/bash

sed -i '/Application/d' out*

for((i=0;i<100;i++))
do
	split -l 5000 out$i.txt out_split$i.txt
done

cat *.txtaa > time_iter.txt
cat *.txtab > time_comp.txt
cat *.txtac > time_comm.txt
