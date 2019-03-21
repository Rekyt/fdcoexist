#!/bin/bash


for i in 1 2 3 4 5 6
do
    echo "Begin $i"
    Rscript --vanilla inst/cluster_bigmem.R $i > big_mem_20190320_$i.out
    wait
    echo "End $i"
done
