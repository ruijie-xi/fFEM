#!/bin/bash

# test prepost

outputdir=output/prepost
mkdir -p $outputdir

smoother=(2) # J, G-S
n_pre=(1 2 5 10 20)
n_post=(1 2 5 10 20)
nlevel=2
n_coarse=64
n_fine=128
for i in {0..4}
do
    outputfile=$outputdir/output_${n_coarse}_${n_fine}_${n_pre[i]}_${n_post[i]}_${smoother}.txt
    ./test_multigrid $nlevel ${n_coarse} ${n_fine} ${n_pre[i]} ${n_post[i]} ${smoother} $outputfile
done

