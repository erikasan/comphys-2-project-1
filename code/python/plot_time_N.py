import numpy as np
from numpy import log10, logspace,log2
import matplotlib.pyplot as plt
import os, sys
import subprocess
import pandas as pd
import seaborn as sns

harr=np.loadtxt("../../output/vmc_time_simpleharmonic.csv",dtype="float",delimiter=",")
number_particles_tot,alphas_tot,energies_vmc,time_vmc=harr.T
harr=np.loadtxt("../../output/numerical_time_simpleharmonic.csv",dtype="float",delimiter=",")
number_particles_tot,alphas_tot,energies_numerical,time_numerical=harr.T
number_particles_tot=(number_particles_tot+1e-10).astype(int)
Ns=[1,5,10,25,50,100]
times=np.zeros(len(Ns)); times_numerical=np.zeros(len(Ns))

for i,N in enumerate(Ns):
    times[i]=np.mean(time_vmc[number_particles_tot==N])
    times_numerical[i]=np.mean(time_numerical[number_particles_tot==N])
print(times_numerical/times)
sns.set_style("darkgrid")
sns.set_context("talk")
plt.plot(Ns,times,"o",label="Analytic 2nd derivative",color="red")
plt.plot(Ns,times_numerical,"o",label="Numerical 2nd derivative",color="blue")

plt.plot(Ns,times,"--",color="red")
plt.plot(Ns,times_numerical,"--",color="blue")

plt.ylabel("Time [ms]")
plt.xlabel("Number of particles")
plt.legend()
plt.xscale("log")
plt.xticks(Ns,Ns)
plt.yscale("log")
plt.tight_layout()
plt.savefig("../../figures/time_numericalvsvmc.pdf")
plt.show()

alphas=np.linspace(0.3,0.7,11)
print("alpha   numerical energy       analytical energy")
for i,alpha in enumerate(alphas):
 print("%5.3f & %10.5f &  %10.5f \\\\hline"%(alpha,energies_numerical[-11+i],energies_vmc[-11+i],))
