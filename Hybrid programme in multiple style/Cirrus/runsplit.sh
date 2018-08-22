#This file is for data processing including cleaning, splitting and gathering the data.
#The data after processing will be stored in time_iter.txt, time_comp.txt and time_comm.txt.
#They store the runtime of each iteration, the computation time in each iteration and the communication time in each iteration respectively.

#!/bin/bash

#sed -i '/time_*/d' out*
#sed -i '/Application/d' out*.txt

for((i=0;i<100;i++))
do
	split -l 5000 out$i.txt out_split$i.txt
done

cat out*.txtaa > time_iter.txt
cat out*.txtab > time_comp.txt
cat out*.txtac > time_comm.txt
