#!/bin/bash


for i in $(seq 1 30)
do
    echo "Begin $i"
    Rscript --vanilla inst/extract_perf.R $i > big_mem_20190326_perf_$i.out
    wait
    echo "End $i"
done
