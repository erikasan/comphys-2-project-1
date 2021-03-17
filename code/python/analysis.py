import os
import sys
# Where to save the figures and data files
DATA_ID = "../../output/"

def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

infile = open(data_path("energies_1616003395561436856.csv"),'r')
import pandas as pd
from pandas import DataFrame

from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, sqrt
from numpy.linalg import inv
import numpy as np
def block(x):
    # preliminaries

    n = len(x)
    d = int(log2(n))
    s, gamma = zeros(d), zeros(d)
    mu = mean(x)

    # estimate the auto-covariance and variances
    # for each blocking transformation
    for i in arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])

    # generate the test observator M_k from the theorem
    M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    return mu, s[k]/2**(d-k)



xy = loadtxt(infile,skiprows=1,usecols=(0),delimiter=",")
maxval=int( log2(len(xy)))
start=2
xvals=np.logspace(start,maxval,maxval-start+1,base=2,dtype=int)
means=np.zeros(len(xvals))
stds=np.zeros(len(xvals))
print(xvals)
for i,xval in enumerate(xvals):
    #print(i,xval)
    #print(xy[0:xval])
    (meany, vary) = block(xy[0:xval])
    std = sqrt(vary)
    means[i]=meany;
    stds[i]=std;
    data ={'Mean':[meany], 'STDev':[std]}
    frame = pd.DataFrame(data,index=['Values'])
    print(frame)
import matplotlib.pyplot as plt
xvals*=100
plt.errorbar(xvals,means,yerr=stds)
plt.xscale("log",base=2);
plt.xlabel("Number of Monte Carlo Cycles")
plt.ylabel("Mean energy")
plt.savefig("test.png")
plt.show()
