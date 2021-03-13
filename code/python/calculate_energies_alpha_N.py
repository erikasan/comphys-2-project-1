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
stepLength=1.0;
equilibration=0.1;
alphas=np.linspace(0.4,0.6,5);#[1/8*np.sqrt(2)**i for i in range(0,10)]
print("%10s %10s %10s %10s"%("Alpha ","Num_part"," Num_dim","E"))
for number_particles in number_particleslist:
    for number_dimensions in number_dimensionslist:
        N = int(1e5)*number_particles
        for alpha in alphas:
            bashCommand="./vmc %d %d %d %f %f %f %f %d"%(number_dimensions,number_particles,N,omega,alpha,stepLength,equilibration,2021)
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd="../cpp/build/",shell=False)
            output, error = process.communicate()
            infile=pd.read_csv(filepath_or_buffer="../../output/sympleharmonic.csv",header=0)
            alpha=np.array(infile["alpha"][-1:])[0]
            energies=np.array(infile["energy"][-1:])[0]
            print("%10.3f %10d %10d %10.3f"%(alpha,number_particles,number_dimensions,energies))
