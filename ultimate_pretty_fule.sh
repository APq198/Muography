#!/bin/bash
numOfThreads=4
prefix="G4W_"
suffix="_$1"
filename="end_output11.txt"
i=0
while [ $i -lt $numOfThreads ]	
do
	##
	##
	##
	input="./$prefix$i$suffix"
	echo $input
	((i++))
done