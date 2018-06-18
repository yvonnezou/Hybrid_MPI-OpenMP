#!/bin/bash

#FILE = $1

for((i=0;i<10;i++))
do
	qsub cfd.pbs
done
