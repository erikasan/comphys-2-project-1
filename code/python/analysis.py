import os
import sys
from blocking import block, naive
import numpy as np
from numpy import log2, sqrt
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
import seaborn as sns
# Where to save the figures and data files
DATA_ID = "../../output/"
def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

filename="energies_100_0.10_0.4830_23_"; num_threads=4


values=np.empty((0));
for i in range(num_threads):
    if num_threads==1:
        infile = open(data_path("%s.csv"%filename),'r')
        values=np.append(values,np.loadtxt(infile,skiprows=1,usecols=(0),delimiter=","))
        infile.close()
        print(len(values)/(i+1))

    else:
        infile = open(data_path("%s%d.csv"%(filename,i)),'r')
        values=np.append(values,np.loadtxt(infile,skiprows=1,usecols=(0),delimiter=","))
        infile.close()
        print(len(values)/(i+1))
        print(values[0:100])
Esquared=values*values
maxval=int( log2(len(values)))
start=2
xvals=np.logspace(start,maxval,maxval-start+1,base=2,dtype=int)
#xvals=[2**maxval]
means=np.zeros(len(xvals))
stds=np.zeros(len(xvals))
naive_standarderror=np.zeros(len(xvals))
#np.random.shuffle(values)
print("MC steps   Mean    standarderror_blocking standarderror_naive")
for i,xval in enumerate(xvals):
    #print(i,xval)
    #print(xy[0:xval])
    (meany, vary) = block(values[0:xval])
    std = sqrt(vary)
    means[i]=meany;
    stds[i]=std;
    naive_standarderror[i]=naive(values[0:xval],Esquared[0:xval])
    print("$2^{%d}$ & %5.6f & %5.6f & %5.5f\\\\ \\hline "%(int(log2(xval)+1e-5),meany,std,naive_standarderror[i]))

sns.set_style("darkgrid")
sns.set_context("talk")
plt.errorbar(xvals,means,yerr=stds,label="Blocking")
plt.errorbar(xvals,means,linestyle="--",yerr=naive_standarderror,label="Naive")
plt.xscale("log",base=2);
plt.xlabel("Number of Monte Carlo Cycles")
plt.ylabel("Mean energy")
plt.legend()
plt.tight_layout()
#plt.savefig("../../figures/blocking.pdf")
#plt.show()
