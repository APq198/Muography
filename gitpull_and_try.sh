#!/bin/bash
git reset --hard HEAD
git pull
#cmake ..
make 
. remove_files.sh output11.txt
###rm end_output11.txt
./main
. ultimate_pretty_file.sh output11.txt
. find_muons.sh end_output11.txt
echo "================"
echo "===== Done ====="
echo "================"