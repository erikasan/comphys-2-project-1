import os
import sys
from blocking import block
import numpy as np
from numpy import log2, sqrt
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt

# Where to save the figures and data files
DATA_ID = "../../output/"

def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

filename="energies_test"; num_threads=4
infile = open(data_path("energies_test.csv"),'r')


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
maxval=int( log2(len(values)))
start=2
xvals=np.logspace(start,maxval,maxval-start+1,base=2,dtype=int)
#xvals=[2**maxval]
means=np.zeros(len(xvals))
stds=np.zeros(len(xvals))
print(xvals)
#np.random.shuffle(values)
for i,xval in enumerate(xvals):
    #print(i,xval)
    #print(xy[0:xval])
    (meany, vary) = block(values[0:xval])
    std = sqrt(vary)
    means[i]=meany;
    stds[i]=std;
    data ={'Mean':[meany], 'STDev':[std]}
    frame = pd.DataFrame(data,index=['Values'])
    print(frame)
plt.errorbar(xvals,means,yerr=stds)
plt.xscale("log",base=2);
plt.xlabel("Number of Monte Carlo Cycles")
plt.ylabel("Mean energy")
plt.savefig("test.png")
plt.show()
