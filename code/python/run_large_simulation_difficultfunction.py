import os
import subprocess
import sys
from blocking import block
import numpy as np
from numpy import log2, sqrt
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
import seaborn as sns
Ns=[3,5,10,25,50,100]
dts=[0.1,0.1,0.1,0.1,0.1,0.1]
alphas=[0.4994,0.4988,0.4975,0.4939,0.489,0.483]
num_threads=4;
number_dimensions=3
numberOfSteps=2**23
seed=1000
for i,number_particles in enumerate(Ns):
    stepLength=dts[i];
    alpha=alphas[i];
    equilibration=1000*number_particles
    filename="%d_%.2f_%.4f_%d_"%(number_particles,stepLength,alpha,int(log2(numberOfSteps)+1e-7))
    bashCommand="./vmc_parallel %d \
    %d %d %f %f %d  %d %s %s\
     %s %d"%(number_dimensions,number_particles,numberOfSteps,alpha,stepLength,equilibration,seed,"EO","IMP",filename,num_threads)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd="../cpp/build/",shell=False)
    output, error = process.communicate()
    print("Finished %d"%number_particles)
