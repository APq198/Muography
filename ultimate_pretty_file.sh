#!/bin/bash
numOfThreads=8
prefix="G4W_"
suffix="_$1"
filename="end_$1"
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
		[[ $var =~ ^u.* ]] && continue
		echo "$var" >> "$filename"
		done < "$input"
	((i++))
done