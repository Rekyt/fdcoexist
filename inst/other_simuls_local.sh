#!/bin/bash

for i in $(seq 1 13)
do
    echo "Begin $i"
    /c/Users/grenie/Documents/R/R-3.5.2/bin/Rscript.exe --arch x64 --vanilla inst/other_simuls.R $i > other_simuls_20190403_$i.out
    wait
    echo "End $i"
done
