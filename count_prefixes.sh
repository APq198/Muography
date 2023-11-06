#!/bin/bash
input="./$1"
counter=0
while IFS= read -r var
do
  #
  # if value of $var starts with "a", ignore it
  #
  ### [[ $var =~ ^a.* ]] && continue
  if [[ $var =~ ^"$2" ]];
  then
	((counter++))
  fi
done < "$input"
echo "Counted $counter lines that start with "$2""