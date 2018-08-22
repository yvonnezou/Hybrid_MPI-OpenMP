#This file is used to upload the t2cfd.pbs file to backend 100 times

#!/bin/bash

for((i=0;i<100;i++))
do
	qsub t2cfd.pbs
done
