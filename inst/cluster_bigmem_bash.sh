#!/bin/bash


for i in $(seq 1 30)
do
    echo "Begin $i"
    Rscript --vanilla inst/cluster_bigmem.R $i > big_mem_20190324_$i.out
    wait
    echo "End $i"
done
