import numpy as np
from numpy import log10, logspace,log2,linspace
import matplotlib.pyplot as plt
import os, sys
import subprocess
import pandas as pd
import seaborn as sns
number_particles=[3,5,10,25,50,100]
number_dimensions=3;
N = int(1e3)
equilibration=int(1e3)
stepLengths=[0.01,0.1,1];
alpha=0.5
num_thread=4
for num_part in number_particles:
    N=1e4*num_part
    equilibration=int(N/10)
    for stepLength in stepLengths:
        bashCommand="./gradientdescent_parallel %d %d %d   %f   %f   %d   %d   %s %s %s %d"%(3,num_part,N,alpha,stepLength,equilibration,1130,"EO","IMP","no",num_thread)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd="../cpp/build/",shell=False)
        output, error = process.communicate()
        print("Done %.2f steplength, %d particles"%(stepLength,num_part))
