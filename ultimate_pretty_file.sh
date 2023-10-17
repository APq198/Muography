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
	while IFS= read -r var
		do
		# if value of $var starts with "a", ignore it
		[[ $var =~ ^a.* ]] && continue
		echo "$var" >> "$filename"
		done < "$input"
	((i++))
done