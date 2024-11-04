#!/bin/bash

#test convergence
outputdir=output/convergence
mkdir -p $outputdir

smoother=2 # G-S
n_pre=5
n_post=5
nlevel=2
n_coarse=(4 8 16 32 64 128)
n_fine=(8 16 32 64 128 256)
for i in {0..5}
do
    outputfile=$outputdir/output_${n_coarse[i]}_${n_fine[i]}_${n_pre}_${n_post}_${smoother}.txt
    ./test_multigrid $nlevel ${n_coarse[i]} ${n_fine[i]} $n_pre $n_post $smoother $outputfile
done