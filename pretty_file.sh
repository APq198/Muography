#!/bin/bash
input="/path/to/txt/file"
while IFS= read -r var
do
  #
  # if value of $var starts with #, ignore it
  #
  [[ $var =~ ^a.* ]] && continue
  echo "$var" >> file2.txt
done < "$input"