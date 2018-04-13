# coding: utf-8

import numpy as np
import params_maooam
from params_maooam import ndim, tw, t_run, t_trans, dt
import integrator
import time
import sys
import aotensor as aotensor_mod
import tl_ad
import scipy
from scipy import linalg, matrix
import matplotlib.pyplot as plt

def print_progress(p):
    sys.stdout.write('Progress {:.2%} \r'.format(p))
    sys.stdout.flush()

class bcolors:
    """to color the instructions in the console"""
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


print (bcolors.OKBLUE + "Starting loading the data file ..." + bcolors.ENDC)

# Time
t = 0.
N = int(t_run/tw) # number of fwd transient steps; t_run is the fwd transient time
test_run = 1.e5 # compute CLVs from t=0 to 1.e5 
test_N = int(test_run/tw) # number of output sets of CLVs from t=0 to 1.e5
btrans_run = 1.e5 # backward transient time
btrans_N = int(btrans_run/tw)
aotensor = aotensor_mod.aotensor

print ("test_N = %d" % test_N)
print ("btrans_N = %d" % btrans_N)

# Initialize variables from data files
q_mat_0 = np.loadtxt('BLV.dat')[1:ndim+1,:] # BLV at Xhist[N+1]
total_run = test_run + btrans_run
total_N = test_N + btrans_N
Xhist = np.loadtxt('evol_field.dat')[(N+1):(N+1)+total_N,1:(ndim+1)] #!!! modify this later
q_mat_old = q_mat_0
t_up = tw/total_run*100
tint = 10. #Default value is the same as tw e.g. 0.1
tint_N = int(np.round(tint/tw))
M = np.identity(ndim)
counter2 = 0
print ("total_N = %d" % total_N)

##### Forward dynamics #####
print (bcolors.OKBLUE + "Starting forward dynamics ..." + bcolors.ENDC)
for counter in np.arange(0,total_N):
    X=Xhist[counter,:]
    M_tmp = tl_ad.compute_tlm(X,tw,aotensor) # compute M[n-1,n] based on X[n-1], n=1,2,...,total_N
    M = np.dot(M_tmp,M)

    if (counter+1) % tint_N == 0:
        M2 = np.dot(M,q_mat_old) 
        q_mat_new, r_mat_new = np.linalg.qr(M2) # generate Q_n, R_n, n=1,2,...,total_N

        # store q_mat_new and r_mat_new in Qhist and Rhist
        if counter2 == 0:
            Qhist = q_mat_new
            Rhist = r_mat_new
            counter2 +=  1
        else:
            Qhist = np.append(Qhist,q_mat_new,axis=0) # Q_n = Q[(n-1)*ndim:n*ndim,:]
            Rhist = np.append(Rhist,r_mat_new,axis=0) # R_n = R[(n-1)*ndim:n*ndim,:]
            counter2 += 1

        q_mat_old = q_mat_new
        M = np.identity(ndim)

    t += tw

    if t/total_run*100 % 0.1 < t_up:
        print_progress(t/total_run)

print (bcolors.OKBLUE + "Complete forward dynamics ..." + bcolors.ENDC)
print (counter2)

# Output Qhist and Rhist as txt data file
np.savetxt('Qhist_CLV.dat',Qhist)
np.savetxt('Rhist_CLV.dat',Rhist)

print (bcolors.OKBLUE + "Start backward transient ..." + bcolors.ENDC)

# Initialize upper triangular matrix C_mat
C_mat = np.zeros([ndim,ndim])
for i in np.arange(0,ndim):
    tmp = np.random.randn(i+1)
    C_mat[0:tmp.shape[0],i] = tmp/np.linalg.norm(tmp)

##### Backward transient #####
for counter in np.arange(0,int(np.round(btrans_N/tint_N))):
    ind = int(np.round(total_N/tint_N)-counter)
    r_mat = Rhist[(ind-1)*ndim:ind*ndim,:] # last R = Rhist[(N-1)*36 : N*36,:]
    r_mat_inv = np.linalg.inv(r_mat)

    cr_tmp = np.dot(r_mat_inv,C_mat) # R_(n+k)^{-1}*C_(n+k)
    # Construct D_mat
    D_mat = np.zeros([ndim,ndim])
    for i in np.arange(0,ndim):
        D_mat[i,i] = 1./np.linalg.norm(cr_tmp[:,i])

    C_mat = np.dot(cr_tmp,D_mat)
    t -= tint

    if t/total_run*100 % 0.1 < t_up:
        print_progress(t/total_run)

print (bcolors.OKBLUE + "Complete backward transient ..." + bcolors.ENDC)
# Note: at this step, ind = total_N-(btrans_N-1), therefore the next step should ind=total_N-btrans_N

##### Backward dynamics #####
counter2 = 0
print (bcolors.OKBLUE + "Start backward dynamics ..." + bcolors.ENDC)
for counter in np.arange(0,int(np.round(test_N/tint_N))):
    ind = int(np.round((total_N-btrans_N)/tint_N)-counter)     
    q_mat = Qhist[(ind-1)*ndim:ind*ndim,:]
    
    if counter2==0:
      CLVhist = np.dot(q_mat,C_mat)
      counter2 += 1
    else:
      CLVhist = np.append(np.dot(q_mat,C_mat),CLVhist,axis=0)
      counter2 += 1

    # Compute the previous C_n = R_(n+k)^{-1} * C_(n+k) * D_mat
    r_mat = Rhist[(ind-1)*ndim:ind*ndim,:]
    r_mat_inv = np.linalg.inv(r_mat)
    cr_tmp = np.dot(r_mat_inv,C_mat) # R_(n+k)^{-1}*C_(n+k)

    D_mat = np.zeros([ndim,ndim])
    for i in np.arange(0,ndim):
        D_mat[i,i] = 1./np.linalg.norm(cr_tmp[:,i])  

    C_mat = np.dot(cr_tmp,D_mat) # C_n
    t -= tint

    if t/total_run*100 % 0.1 < t_up:
        print_progress(t/total_run)

# Don't forget to save CLV at t=0
CLVhist = np.append(np.dot(q_mat_0,C_mat),CLVhist,axis=0)
print (bcolors.OKBLUE + "Complete backward dynamics ..." + bcolors.ENDC)

# Save CLVhist
np.savetxt('CLVhist.dat',CLVhist)

# Save the 1st C matrix
np.savetxt('1st_Cmat.dat',C_mat)

