#!/bin/bash

# test smoother

outputdir=output/coarse
mkdir -p $outputdir

smoother=2 # G-S
n_pre=5
n_post=5
nlevel=2
n_coarse=(30 40 50 60 70)
n_fine=100
for i in {0..4}
do
    outputfile=$outputdir/output_${n_coarse[i]}_${n_fine}_${n_pre}_${n_post}_${smoother}.txt
    ./test_multigrid $nlevel ${n_coarse[i]} ${n_fine} $n_pre $n_post ${smoother} $outputfile
done

