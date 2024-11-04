#!/bin/bash

# test smoother

outputdir=output/smoother
mkdir -p $outputdir

smoother=(1 2) # J, G-S
n_pre=5
n_post=5
nlevel=2
n_coarse=64
n_fine=128
for i in {0..1}
do
    outputfile=$outputdir/output_${n_coarse}_${n_fine}_${n_pre}_${n_post}_${smoother[i]}.txt
    ./test_multigrid $nlevel ${n_coarse} ${n_fine} $n_pre $n_post ${smoother[i]} $outputfile
done

