import numpy as np
import params_maooam
from params_maooam import ndim, natm, noc, tw, t_run, t_trans, dt, f0
import integrator
import time
import sys
import aotensor as aotensor_mod
import tl_ad_solo

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

print (bcolors.OKBLUE + "Starting the time evolution ..." + bcolors.ENDC)
t = 0.
t_up = dt/t_run*100
N = int(np.round(t_run/tw))
solo = 'ocn'
Xhist = np.loadtxt('evol_field_solo.dat')[0:(N+1),1:ndim+1]
T = time.clock()

X=Xhist[0,:]

aotensor = aotensor_mod.aotensor
I = np.identity(ndim)
counter = 0
err = 1.

X_atm = X[0:2*natm]
X_ocn = X[2*natm:ndim] 
M_solo = tl_ad_solo.compute_tlm_solo(X_atm,X_ocn,tw,aotensor,solo) 
LE_ave = np.zeros([1,2*noc])
q_mat_new,r_mat_new = np.linalg.qr(M_solo) 
LE_sum = np.log(np.abs(r_mat_new.diagonal()))/tw
r_mat_prod = r_mat_new

#while t < t_run:
while counter < N:

    t += tw  
    X = Xhist[counter+1,:]
    X_atm = X[0:2*natm]
    X_ocn = X[2*natm:ndim] 
    
    M_solo = tl_ad_solo.compute_tlm_solo(X_atm,X_ocn,tw,aotensor,solo)    
    M2 = np.dot(M_solo,q_mat_new)
    q_mat_new,r_mat_new = np.linalg.qr(M2)    
    LE_sum = LE_sum + np.log(np.abs(r_mat_new.diagonal()))/tw
    r_mat_prod = np.dot(r_mat_new,r_mat_prod)
       
    if counter == 0:
        LE_tmp = LE_sum/(counter+2)
        sort_inds = LE_tmp.argsort()
        LE_ave[counter,:] = LE_tmp[sort_inds[::-1]]
    else:
        LE_tmp = LE_sum/(counter+2)
        sort_inds = LE_tmp.argsort()
        LE_ave = np.append(LE_ave,[LE_tmp[sort_inds[::-1]]],axis=0)
        
    counter +=1
    if t/t_run*100 % 0.1 < t_up:
        print_progress(t/t_run)

# conert LEs from unit model time to unit day
LE_ave = LE_ave*(86400*f0)
LE_unsort = LE_sum/(counter+1)*8.9
BLV = q_mat_new 

print (bcolors.OKBLUE + "Evolution finished " + bcolors.ENDC)
print (bcolors.OKBLUE + "Time clock :" + bcolors.ENDC)
print (time.clock()-T)

# save the BLVs as the row vectors in BLV.dat file
fichier = open("BLV_4_104_%s.dat" % solo, "w")
for i in np.arange(0,2*noc):
    fichier.write(str(LE_unsort[i])+" ")
fichier.write("\n")
for j in np.arange(0,2*noc):
    for i in np.arange(0,2*noc):
        fichier.write(str(BLV[j,i])+" ")
    fichier.write("\n")

fichier.close()

print (np.shape(LE_ave))
print (counter)
print (LE_ave[-1,:])
#print (np.linalg.norm(LE_ave[-1,1:ndim+1]-LE_ave[-2,1:ndim+1]))
