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
Xhist = np.loadtxt('../fortran_da/nature.dat')[0:(N+1),1:ndim+1]

ens_num = sys.argv[1]
tw_da_s = sys.argv[2]
tw_da = float(tw_da_s)

infl = 1.0
#gainfile = '../fortran_da/gain_etkf_%s_%3.1f%6.1E.dat' % (str(ens_num),infl,tw_da)
gainfile = '../fortran_da/gain_etkf_%s_%3.1f%s.dat' % (str(ens_num),infl,tw_da_s)
Khist = np.loadtxt(gainfile) # the same as the number of DA cycles
Hmat = np.identity(2*natm) # Assume perfect observations
print("Khist.shape = ", Khist.shape)
solo = 'atm'
T = time.clock()

X=Xhist[0,:]

aotensor = aotensor_mod.aotensor
I = np.identity(ndim)
counter = 0
err = 1.

X_atm = X[0:2*natm]
X_ocn = X[2*natm:ndim] 

LE_ave = np.zeros(2*natm)
LE_sum = np.zeros(2*natm)
LEhist = np.zeros([1,2*natm])
q_mat_new = np.identity(2*natm)
r_mat_new = np.identity(2*natm)
LEhist[0,:] = LE_sum
r_mat_prod = np.identity(2*natm)
cnt_da = 0

#while t < t_run:
while counter < N:

    t += tw  
    X = Xhist[counter,:]
    X_atm = X[0:2*natm]
    X_ocn = X[2*natm:ndim] 
    
    if counter == 0: 
      M_solo = tl_ad_solo.compute_tlm_solo(X_atm,X_ocn,tw+dt,aotensor,solo)
    else:
      M_solo = tl_ad_solo.compute_tlm_solo(X_atm,X_ocn,tw,aotensor,solo)

    if np.mod(counter+1,np.int(tw_da/tw)) == 0 :
        Kmat = Khist[2*natm*(cnt_da):2*natm*(cnt_da+1),:] 
        M_da = np.dot(np.identity(2*natm)-np.dot(Kmat,Hmat),M_solo)   
        cnt_da += 1
    else:
        M_da = M_solo
   
    M2 = np.dot(M_da,q_mat_new)
    q_mat_new,r_mat_new = np.linalg.qr(M2)    
    LE_sum = LE_sum + np.log(np.abs(r_mat_new.diagonal()))/tw
    r_mat_prod = np.dot(r_mat_new,r_mat_prod)
       
    
    LE_tmp = LE_sum/(counter+1)
    sort_inds = LE_tmp.argsort()
    LE_ave = LE_tmp[sort_inds[::-1]]
    
    if counter == 0:
        LEhist[counter,:] = LE_ave
    else:
        LEhist = np.append(LEhist,[LE_ave],axis=0)    

    counter +=1
    if t/t_run*100 % 0.1 < t_up:
        print_progress(t/t_run)

# conert LEs from unit model time to unit day
LE_ave = LE_ave*(86400*f0)
LE_unsort = LE_sum/(counter+1)*(86400*f0)
BLV = q_mat_new 

print (bcolors.OKBLUE + "Evolution finished " + bcolors.ENDC)
print (bcolors.OKBLUE + "Time clock :" + bcolors.ENDC)
print (time.clock()-T)

# save LEhist
np.savetxt('LEhist_etkf_%s_%3.1f%s.dat' % (str(ens_num),infl,tw_da),LEhist*(86400*f0))

# save the BLVs as the row vectors in BLV.dat file
fichier = open("BLV_5_104_%s_etkf_%s%s.dat" % (solo,str(ens_num),tw_da_s), "w")
for i in np.arange(0,2*natm):
    fichier.write(str(LE_ave[i])+" ")
fichier.write("\n")
for j in np.arange(0,2*natm):
    for i in np.arange(0,2*natm):
        fichier.write(str(BLV[j,i])+" ")
    fichier.write("\n")

fichier.close()

print (np.shape(LE_ave))
print (counter)
#print (exp)
#print (LE_ave[-1,:])
#print (np.linalg.norm(LE_ave[-1,1:ndim+1]-LE_ave[-2,1:ndim+1]))
