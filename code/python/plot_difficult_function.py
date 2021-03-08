import numpy as np
from numpy import log10, logspace,log2
import matplotlib.pyplot as plt
import os, sys
import subprocess
import pandas as pd
import seaborn as sns

number_particles=10
number_dimensions=3;
N = int(1e5)
omega=1
stepLength=0.1;
equilibration=0.1;
beta=2.82843
a=0.0043
seed=12;
alphas=[1/8*np.sqrt(2)**i for i in range(0,10)]
for alpha in alphas:
    bashCommand="./monster %d %d %d   %f %f  %f  %f  %f       %f     %d"%(number_dimensions,number_particles,N,omega,alpha,stepLength,equilibration,a,beta, seed)
    #bashCommand="./monster 3 10 100000 1 0.7 0.1 0.1 0.0043 2.82843 12345"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd="../cpp/build/",shell=False)
    output, error = process.communicate()
    print(output)
infile=pd.read_csv(filepath_or_buffer="../../output/sympleharmonic.csv",header=0)
print(infile)
alphas=infile["alpha"][-len(alphas):]
energies=infile["energy"][-len(alphas):]
kinetic_energies=infile["kin_en"][-len(alphas):]
potential_energies=infile["pot_en"][-len(alphas):]

sns.set()
plt.plot(alphas,energies,label=r"$E_{TOT}$")
plt.plot(alphas,potential_energies,label=r"$V$")
plt.plot(alphas,kinetic_energies,label=r"$K$")
plt.xlabel(r"$\alpha$")
plt.ylabel(r"E (Hartree)")
plt.legend()
plt.grid()
plt.show()
