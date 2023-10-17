#!/bin/bash
numOfThreads=4
prefix="G4W_"
suffix="_$1"
filename="end_output11.txt"
for i in {0..$numOfThreads-1}
do
	##
	##
	##
	input="./"+prefix+$i+suffix
	echo input
done