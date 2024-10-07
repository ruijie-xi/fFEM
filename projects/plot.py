import matplotlib.pyplot as plt
import numpy as np

N = 100
coeff_list = [10,20,30,40,50]
smooth_list = [1, 2]

cut_step = 100

for smooth in smooth_list:
    
    plt.figure()
    
    smooth_name = "Jacobi" if smooth == 1 else "Gauss-Seidel"
    
    plt.title("Residue vs. Iteration step ({0} smoothing)".format(smooth_name))

    for coeff in coeff_list:
        inputdir = "N{0}_coeff{1}_smooth{2}".format(N, coeff, smooth)
        with open(inputdir + "/residue.dat") as f:
            lines = f.readlines()
            step = []
            res = []
            res0 = float(lines[0].split()[1])
            for line in lines[:cut_step+1]:
                step.append(int(line.split()[0]))
                res.append(float(line.split()[1])/res0)
            plt.plot(step, res, label="k={0}".format(coeff))
            
    plt.legend()
    plt.xlabel("Iteration step")
    plt.ylabel("Residue")
    
    plt.yscale("log")
    
    plt.savefig(smooth_name + "_residue.png")
    
    plt.close()
    
    
    
    
    