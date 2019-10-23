# This code is used for rmse comparison in different experiments

import numpy as np
import matplotlib.pyplot as plt

true = np.loadtxt('nature.dat')[1:1000001,1:39]
free = np.loadtxt('freerun.dat')[1:1000001,1:39]
etkf = np.loadtxt('Xam_etkf_atm_37_1.0.2E+01_01.dat')[0:1000000,1:39]
etkf2 = np.loadtxt('../fortran_da/Xam_etkf_atm_37_1.0.2E+01_01.dat')[0:1000000,1:37]
etkf3 = np.loadtxt('../fortran_da/Xam_etkf_cpl_37_1.0.2E+01_01.dat')[0:1000000,1:37]

tmf_atm  = np.linalg.norm(true[:,0:20]-free[:,0:20],axis=1)
tma_atm  = np.linalg.norm(true[:,0:20]-etkf[:,0:20],axis=1)
tma_atm2 = np.linalg.norm(true[:,0:20]-etkf2[:,0:20],axis=1)
tma_atm3 = np.linalg.norm(true[:,0:20]-etkf3[:,0:20],axis=1)

tmf_ocn  = np.linalg.norm(true[:,20:36]-free[:,20:36],axis=1)
tma_ocn  = np.linalg.norm(true[:,20:36]-etkf[:,20:36],axis=1)
tma_ocn2 = np.linalg.norm(true[:,20:36]-etkf2[:,20:36],axis=1)
tma_ocn3 = np.linalg.norm(true[:,20:36]-etkf3[:,20:36],axis=1)

ntime = true.shape[0]
time  = np.arange(1,ntime+1)/8.9

plt.figure(1)
ax1=plt.subplot(2,1,1)
plt.plot(time,tmf_atm,label="free error")
plt.plot(time,tma_atm,label="obs:atm+drf")
plt.plot(time,tma_atm2,label="obs:atm")
plt.plot(time,tma_atm3,label="obs:atm+ocn")
plt.xscale('log')
plt.yscale('log')
#plt.xlabel('time (days)')
plt.ylabel('rmse')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
plt.title('atm states')
#plt.legend()

ax2=plt.subplot(2,1,2)
plt.plot(time,tmf_ocn,label="free error")
plt.plot(time,tma_ocn,label="obs:atm+drf")
plt.plot(time,tma_ocn2,label="obs:atm")
plt.plot(time,tma_ocn3,label="obs:atm+ocn")
plt.xscale('log')
plt.yscale('log')
plt.xlabel('time (days)')
plt.ylabel('rmse')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
plt.title('ocn states')
art2=[]
lgd2=plt.legend(loc='upper center',bbox_to_anchor=(0.5,-0.15),ncol=4,frameon=False)
art2.append(lgd2)
plt.tight_layout()

plt.savefig('rmse_atm_compare.eps')
plt.show()
