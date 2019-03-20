echo "Begin all"
Rscript --vanilla inst/cluster_bigmem.R 1 > big_mem_20190320_1.out
wait
echo "End first; begin second"
Rscript --vanilla inst/cluster_bigmem.R 2 > big_mem_20190320_2.out
wait
echo "End second; begin third"
Rscript --vanilla inst/cluster_bigmem.R 3 > big_mem_20190320_3.out
echo "End all"
