import numpy as np
from numpy import log10, logspace,log2,linspace
import matplotlib.pyplot as plt
import os, sys
import subprocess
import pandas as pd
import seaborn as sns

number_particles=[1,5,10,25,50,100]
number_dimensions=3;
N = int(1e4)
omega=1
stepLength=1;
alphas=linspace(0.3,0.7,11)
totnum=len(alphas)*len(number_particles)
number_particles_tot=np.repeat(number_particles,len(alphas),axis=0)
alphas_tot=np.tile(alphas,len(number_particles))
for num_part in number_particles:
    number_runs=N*num_part
    equilibration=int(0.1*number_runs)
    for alpha in alphas:
        bashCommand="./vmc %d %d %d %f %f %d  %d %s %s %s"%(number_dimensions,num_part,number_runs,alpha,stepLength,equilibration,2021,"HO","VMC","no")
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd="../cpp/build/",shell=False)
        output, error = process.communicate()
        print("Done alpha=%f num_part=%d"%(alpha,num_part))
infile=pd.read_csv(filepath_or_buffer="../../output/sympleharmonic.csv",header=0)
totnum=len(alphas)*len(number_particles)
alphas_vmc=np.array(infile["alpha"][-totnum:])
energies_vmc=np.array(infile["energy"][-totnum:])
kinetic_energies_vmc=np.array(infile["kin_en"][-totnum:])
potential_energies_vmc=np.array(infile["pot_en"][-totnum:])
time_vmc=np.array(infile["time"][-totnum:])



a = np.asarray([number_particles_tot,alphas_tot,energies_vmc,time_vmc]).T
np.savetxt("../../output/vmc_time_simpleharmonic.csv", a, delimiter=",")

for num_part in number_particles:
    number_runs=N*num_part
    equilibration=int(0.1*number_runs)
    for alpha in alphas:
        bashCommand="./vmc %d %d %d %f %f %d  %d %s %s %s"%(number_dimensions,num_part,number_runs,alpha,stepLength,equilibration,2021,"NHO","VMC","no")
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd="../cpp/build/",shell=False)
        output, error = process.communicate()
        print("Done alpha=%f num_part=%d"%(alpha,num_part))
infile=pd.read_csv(filepath_or_buffer="../../output/sympleharmonic.csv",header=0)
totnum=len(alphas)*len(number_particles)
alphas_vmc=np.array(infile["alpha"][-totnum:])
energies_vmc=np.array(infile["energy"][-totnum:])
kinetic_energies_vmc=np.array(infile["kin_en"][-totnum:])
potential_energies_vmc=np.array(infile["pot_en"][-totnum:])
time_vmc=np.array(infile["time"][-totnum:])
a = np.asarray([number_particles_tot,alphas_tot,energies_vmc,time_vmc]).T
np.savetxt("../../output/numerical_time_simpleharmonic.csv", a, delimiter=",")
