#!/bin/bash

outputdir=output/test
mkdir -p $outputdir

smoother=3 # 1 for Jacobi, 2 for Gauss-Seidel
n_pre=5
n_post=5
nlevel=2
n_coarse=50
n_fine=100

sigma=1.0
nu=0.1
beta=1.0

outputfile=$outputdir/output_${n_coarse}_${n_fine}_Pr${n_pre}_Po${n_post}_S${smoother}_sig${sigma}_nu${nu}_beta${beta}.txt

./test_multigrid $nlevel ${n_coarse} ${n_fine} $n_pre $n_post ${smoother} $outputfile $sigma $nu $beta


