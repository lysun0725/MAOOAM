import numpy as np
import params_maooam
from params_maooam import ndim, natm, noc, tw, t_trans, dt, f0
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
t_run = 9.e4
t_up = dt/t_run*100
N = int(np.round(t_run/tw))
Xhist = np.loadtxt('../fortran_da/nature.dat')[0:(N+1),1:ndim+1]

ens_num = sys.argv[1]
tw_da_s = sys.argv[2]
tw_da = float(tw_da_s)
tw_da = 2.5
tw_solo = sys.argv[3]

infl = 1.0
#gainfile = '../fortran_da/gain_etkf_%s_%3.1f%6.1E.dat' % (str(ens_num),infl,tw_da)
gainfile = '../../MAOOAM_solo_atm/fortran_da/h%s/gain_etkf_%s_%3.1f%s.dat' % (tw_solo, str(ens_num),infl,tw_da_s)
gainfile2 = '../../MAOOAM_solo_ocn/fortran_da/h%s/gain_etkf_%s_%3.1f%s.dat' % (tw_solo, str(ens_num),infl,tw_da_s)
Khist_a = np.loadtxt(gainfile) # the same as the number of DA cycles
Khist_o = np.loadtxt(gainfile2)

Hmat = np.identity(ndim)

print("Khist_a.shape = ", Khist_a.shape)
print("Khist_o.shape = ", Khist_o.shape)

T = time.clock()

X=Xhist[0,:]

aotensor = aotensor_mod.aotensor
I = np.identity(ndim)
counter = 0
err = 1.

LE_ave = np.zeros(ndim)
LE_sum = np.zeros(ndim)
LEhist = np.zeros([1,ndim])
Xatm = np.zeros(2*natm)
Xocn = np.zeros(2*noc)
M_solo = np.zeros([ndim,ndim])
Kmat = np.zeros([ndim,ndim])
q_mat_new = np.identity(ndim)
r_mat_new = np.identity(ndim)
LEhist[0,:] = LE_sum
r_mat_prod = np.identity(ndim)
cnt_da = 0

#while t < t_run:
while counter < N:

    t += tw  
    X = Xhist[counter,:]
    Xatm = X[0:2*natm]
    Xocn = X[2*natm:ndim]
    
    if counter == 0: 
      M_solo[0:2*natm,0:2*natm] = tl_ad_solo.compute_tlm_solo(Xatm,Xocn,tw+dt,aotensor,'atm')
      M_solo[2*natm:ndim,2*natm:ndim] = tl_ad_solo.compute_tlm_solo(Xatm,Xocn,tw+dt,aotensor,'ocn')
    else:
      M_solo[0:2*natm,0:2*natm] = tl_ad_solo.compute_tlm_solo(Xatm,Xocn,tw,aotensor,'atm')
      M_solo[2*natm:ndim,2*natm:ndim] = tl_ad_solo.compute_tlm_solo(Xatm,Xocn,tw,aotensor,'ocn')

    if np.mod(counter+1,np.int(tw_da/tw)) == 0 :
        Kmat[0:2*natm,0:2*natm] = Khist_a[2*natm*(cnt_da):2*natm*(cnt_da+1),:]
        Kmat[2*natm:ndim,2*natm:ndim] =  Khist_o[2*noc*(cnt_da):2*noc*(cnt_da+1),:]
        M_da = np.dot(np.identity(ndim)-np.dot(Kmat,Hmat),M_solo)   
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
np.savetxt('LEhist_etkf_%s_%s_%3.1f%s.dat' % ('ucp',str(ens_num),infl,tw_da),LEhist*(86400*f0))

# save the BLVs as the row vectors in BLV.dat file
fichier = open("BLV_5_104_%s_etkf_%s%s.dat" % ('ucp',str(ens_num),tw_da_s), "w")
for i in np.arange(0,ndim):
    fichier.write(str(LE_ave[i])+" ")
fichier.write("\n")
for j in np.arange(0,ndim):
    for i in np.arange(0,ndim):
        fichier.write(str(BLV[j,i])+" ")
    fichier.write("\n")

fichier.close()

print (np.shape(LE_ave))
print (counter)
#print (exp)
#print (LE_ave[-1,:])
#print (np.linalg.norm(LE_ave[-1,1:ndim+1]-LE_ave[-2,1:ndim+1]))
