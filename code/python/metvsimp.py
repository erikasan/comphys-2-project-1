import numpy as np
from numpy import log10, logspace,log2
import matplotlib.pyplot as plt
import os, sys
import subprocess
import pandas as pd
import seaborn as sns

number_particles=10
number_dimensions=3;
equilibration=0.1;
omega=1
Ns=np.array([1e2,1e3,1e4,1e5,1e6,1e7])
alpha=1
steps=[0.01,0.1,0.5,0.7,1]
seeds=np.array(np.linspace(1,100,10),dtype=int)
anal_en=(alpha+(1-4*alpha**2)/(8*alpha))*number_particles*number_dimensions
energies=np.zeros((len(Ns),2*len(steps)),dtype=float)
for j, N in enumerate(Ns):
    for i, stepLength in enumerate(steps):
        equilibration=int(0.1*N)
        for seed in seeds:
            bashCommand="./vmc %d %d %d %f %f %d %d %s %s %s"%(number_dimensions,number_particles,N,alpha,stepLength,equilibration,seed,"HO","VMC","no")
            print(bashCommand)
            print("N: %d steplength: %.2f"%(N, stepLength))
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd="../cpp/build/",shell=False)
            output, error = process.communicate()
            infile=pd.read_csv(filepath_or_buffer="../../output/sympleharmonic.csv",header=0)
            energy=np.array(infile["energy"])[-1]
            energies[j,i]+=energy;
        print("Steplength: %.2f, N: %d, energy: %.3f"%(stepLength,N,energies[j,i]/len(seeds)))
for j, N in enumerate(Ns):
    for i, stepLength in enumerate(steps):
        equilibration=int(0.1*N)
        for seed in seeds:
            bashCommand="./vmc %d %d %d %f %f %d  %d %s %s %s"%(number_dimensions,number_particles,N,alpha,stepLength,equilibration,seed,"HO","IMP","no")
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd="../cpp/build/",shell=False)
            output, error = process.communicate()
            infile=pd.read_csv(filepath_or_buffer="../../output/sympleharmonic.csv",header=0)
            energy=np.array(infile["energy"])[-1]
            energies[j,len(steps)+i]+=energy;
        print("Steplength: %.2f, N: %d, energy: %.3f"%(stepLength,N,energies[j,len(steps)+i]/len(seeds)))
energies/=len(seeds)
energies-=anal_en
energies/=anal_en
energies=np.abs(energies)
sns.set_style("darkgrid")
sns.set_context("talk")
fig, (ax1, ax2) = plt.subplots(1, 2,sharex=True, sharey=True,figsize=(15,8))
for i, stepLength in enumerate(steps):
    ax1.plot(Ns,energies[:,i],label=r"$\Delta x$=%.2f"%(stepLength))
for i, stepLength in enumerate(steps):
    ax2.plot(Ns,energies[:,i+len(steps)],label=r"$\Delta t$=%.2f"%(stepLength))


ax1.legend()
ax2.legend()
"""
ax1.set_ylabel("Relative error")
ax2.set_ylabel("Relative error")

ax1.set_xlabel("Number of simulation steps")
ax2.set_xlabel("Number of simulation steps")
"""
ax1.set_title("Brute Force")
ax2.set_title("Importance Sampling")

ax1.set_xscale("log")
ax1.set_yscale("log")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax1.set_ylim((1e-5,1e-0))
ax2.set_ylim((1e-5,1e-0))
#plt.tight_layout()
fig.text(0.5, 0.01, 'Number of simulation steps', ha='center')
fig.text(0.04, 0.5, 'Relative error', va='center', rotation='vertical')
plt.savefig("../../figures/comparison.pdf")
plt.show()
