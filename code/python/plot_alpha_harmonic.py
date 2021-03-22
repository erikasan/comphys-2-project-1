import numpy as np
from numpy import log10, logspace,log2
import matplotlib.pyplot as plt
import os, sys
import subprocess
import pandas as pd
import seaborn as sns

number_particles=10
number_dimensions=3;
N = int(1e4)*number_particles
omega=1
stepLength=1.0;
equilibration=int(0.1*N);
alphas=np.linspace(0.1,10,10);#[1/8*np.sqrt(2)**i for i in range(0,10)]
for alpha in alphas:
    bashCommand="./vmc %d %d %d %f %f %d  %d %s %s %s"%(number_dimensions,number_particles,N,alpha,stepLength,equilibration,2021,"HO","VMC","test")
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd="../cpp/build/",shell=False)
    output, error = process.communicate()
infile=pd.read_csv(filepath_or_buffer="../../output/sympleharmonic.csv",header=0)
print(infile)
alphas=infile["alpha"][-len(alphas):]
energies=infile["energy"][-len(alphas):]
kinetic_energies=infile["kin_en"][-len(alphas):]
potential_energies=infile["pot_en"][-len(alphas):]

sns.set_style("darkgrid")
sns.set_context("talk")
plt.plot(alphas,energies,label=r"$\langle E_L\rangle$")
plt.plot(alphas,potential_energies,label=r"$V$")
plt.plot(alphas,kinetic_energies,label=r"$K$")
plt.plot(np.array(alphas)[np.argmin(np.array(energies))],np.min(np.array(energies)),"o",color="red")
plt.xlabel(r"$\alpha$")
#plt.xscale("log")
plt.ylabel(r"$E$ ($\hbar\omega$)")
plt.legend()
plt.tight_layout()
#plt.savefig("../../figures/energy_alpha.pdf")
plt.show()
