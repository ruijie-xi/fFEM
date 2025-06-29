#!/bin/bash

N=100
coeff_list=(10 20 30 40 50)
maxsteps=100
smooth_list=(1 2)

for coeff in ${coeff_list[@]}
do
    for smooth in ${smooth_list[@]}
    do
        outputdir="output/N${N}_coeff${coeff}_smooth${smooth}"
        mkdir -p $outputdir
        echo "Running test_jacobi_GS with N=$N, coeff=$coeff, maxsteps=$maxsteps, smoothtype=$smooth"
        ./test_jacobi_GS $N $N $coeff $maxsteps $outputdir $smooth
    done
done

python3 ./plot.py