import numpy as np
from numpy import log10, logspace,log2, pi,exp,sqrt
import matplotlib.pyplot as plt
import os, sys
import subprocess
import pandas as pd
import seaborn as sns
def fx(x):
     return exp(-x*x)/sqrt(pi) #Normalized
def rhof(x):
 return 4*x*x*exp(-x*x)/sqrt(pi) # Normalized


inverse_steplength=20
filename="alpha_0"
infile=np.loadtxt("../../output/positions_%s.csv"%filename,dtype="float",delimiter=",",skiprows=1)
x=infile[:,0]
rho=infile[:,2]
rho_density=infile[:,3]
x_density=infile[:,1]

filename="alpha_0483"
infile=np.loadtxt("../../output/positions_%s.csv"%filename,dtype="float",delimiter=",",skiprows=1)
x_a0043=infile[:,0]
rho_a0043=infile[:,2]
rho_density_a0043=infile[:,3]
x_density_a0043=infile[:,1]

x_density=x_density/np.sum(x_density)*inverse_steplength
rho_density=rho_density/np.sum(rho_density)*inverse_steplength

x_density_a0043=x_density_a0043/np.sum(x_density_a0043)*inverse_steplength
rho_density_a0043=rho_density_a0043/np.sum(rho_density_a0043)*inverse_steplength

sns.set_style("darkgrid")
sns.set_context("talk")
fig, (ax1, ax2) = plt.subplots(1, 2,sharex=False, sharey=True,figsize=(15,8))
"""
ax1.plot(rho,rhof(rho),label=r"$\rho(r)$ analytical",linewidth=3.0)
ax1.plot(rho,rho_density,linestyle="dashed",label=r"$\rho(r)$ numerical",color="red")
"""
ax1.plot(rho,rhof(rho),label=r"$\rho(r)$ for a SHO")
ax1.plot(rho,rho_density,linestyle="dashed",label=r"$\rho(r)$ ; $a=0$",color="red")
ax1.plot(rho_a0043,rho_density_a0043,linestyle="dashed",label=r"$\rho(r)$ ; $a=0.0043$",color="black")

ax1.legend()
ax1.set_xlabel(r"$r/a_0$")

ax1.set_title(r"$\rho(r)$")
"""
ax2.plot(x,fx(x),label=r"$f_x(x)$ analytical",linewidth=3.0)
ax2.plot(x,x_density,linestyle="dashed",label=r"$f_x(x)$ numerical",color="red")
"""
ax2.plot(x,fx(x),label=r"$f_x(x)$ for a SHO")
ax2.plot(x,x_density,linestyle="dashed",label=r"$f_x(x)$ ; $a=0$",color="red")
ax2.plot(x_a0043,x_density_a0043,linestyle="dashed",label=r"$f_x(x)$ ; $a=0.0043$",color="black")

ax2.legend()
ax2.set_xlabel(r"$x/a_0$")

ax2.set_title(r"$f_x(x)$")

ax1.set_xlim([0,3.5])
ax2.set_xlim([-2.5,2.5])

fig.text(0.04, 0.5, 'probability/density', va='center', rotation='vertical')
plt.savefig("../../figures/distributionfunction.pdf")
plt.show()
