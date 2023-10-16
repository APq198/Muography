#!/bin/bash
input="./test.txt"
while IFS= read -r line
do
  echo "$line"
done < "$input"