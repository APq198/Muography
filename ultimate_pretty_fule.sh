#!/bin/bash
numOfThreads=4
prefix="G4W_"
suffix="_$1"
filename="end_output11.txt"
for i in {1..$numOfThreads}
do
	##
	##
	##
	echo $i
	let j=$i-1
	input="./$prefix$i$suffix"
	echo $input
done