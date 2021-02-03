import numpy as np
import matplotlib.pyplot as plt
import os, sys
import subprocess
number_particles=10
number_dimensions=3;
N = int(1e6)
omega=10
stepLength=0.1;
equilibration=0.1;
alphas=np.linspace(1/8*omega,2*omega,10)
for alpha in alphas:
    bashCommand="./vmc %d %d %d %f %f %f %f"%(number_dimensions,number_particles,N,omega,alpha,stepLength,equilibration)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd="../cpp/build/",shell=False)
    output, error = process.communicate()
infile=np.loadtxt("../../output/sympleharmonic.csv",skiprows=1, dtype="float",delimiter=",",usecols=range(1,11))
