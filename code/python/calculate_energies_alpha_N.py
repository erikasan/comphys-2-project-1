import numpy as np
from numpy import log10, logspace,log2,linspace
import matplotlib.pyplot as plt
import os, sys
import subprocess
import pandas as pd
import seaborn as sns

number_particleslist=[1,10,100,500]
number_dimensionslist=[1,2,3];

omega=1
stepLength=0.7;
equilibration=10000;
alphas=np.linspace(0.4,0.6,5);#[1/8*np.sqrt(2)**i for i in range(0,10)]
for alpha in alphas:
    print(r"& $\alpha=%.2f$"%alpha, end=" ")
print ("\\\\ \\hline")
for number_particles in number_particleslist:
    for number_dimensions in number_dimensionslist:
        N = int(1e5)
        print(r"$\langle E_L \rangle$, N=%d, %dD &"%(number_particles,number_dimensions),end="")
        for alpha in alphas:
            bashCommand="./vmc %d %d %d %f %f %d  %d %s %s %s"%(number_dimensions,number_particles,N,alpha,stepLength,equilibration,2021,"HO","IMP","no")
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd="../cpp/build/",shell=False)
            output, error = process.communicate()
            infile=pd.read_csv(filepath_or_buffer="../../output/sympleharmonic.csv",header=0)
            alpha=np.array(infile["alpha"][-1:])[0]
            energies=np.array(infile["energy"][-1:])[0]
            if alpha==alphas[-1]:
                print("%10.3f"%energies,end=" ")
            else:
                print("%10.3f &"%energies,end=" ")
        print("\\\\ \\hline")
