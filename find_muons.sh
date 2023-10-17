#!/bin/bash
input="./$1"
counter=0
while IFS= read -r var
do
  #
  # if value of $var starts with "a", ignore it
  #
  ### [[ $var =~ ^a.* ]] && continue
  if [[$var =~ ^"mu"]];
  then
	((counter++))
  fi
done < "$input"
echo "Counted: $counter muons ( prefix "mu" )"