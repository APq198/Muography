#!/bin/bash
input="./$1"
while IFS= read -r var
do
  #
  # if value of $var starts with "a", ignore it
  #
  [[ $var =~ ^a.* ]] && continue
  echo "$var" >> NEW_$1
done < "$input"