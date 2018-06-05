import numpy as np
import matplotlib.pyplot as plt
import sys

ens_num = sys.argv[1]
tw_da = sys.argv[2]
tw_solo = sys.argv[3]
infl = 1.0
da_file = "Xam_etkf_ocn_%s_%3.1f%s.dat" % (str(ens_num),infl,tw_da)
#sp_file = "sprd_etkf_ocn_%s_%3.1f%s.dat" % (str(ens_num),infl,tw_da)
sp_file = "../../MAOOAM_solo_ocn/fortran_da/h%s/Xam_etkf_%s_%3.1f%s.dat" % (tw_solo,str(ens_num),infl,tw_da)

true = np.loadtxt('nature.dat')[1:900002,21:37]
#free = np.loadtxt('freerun_atm.010.dat')[1:900002,1:37]
etkf = np.loadtxt(da_file)[:,21:37]
sprd = np.loadtxt(sp_file)[:,1:17]
days = np.arange(1,900002)*0.1/8.9

#err1 = np.linalg.norm(true-free,axis=1)/np.sqrt(36)
err2 = np.linalg.norm(true-etkf,axis=1)/np.sqrt(16)
err3 = np.linalg.norm(true-sprd,axis=1)/np.sqrt(16)

#plt.plot(err1,label = 'true-free')
plt.plot(days,err2,label = 'true-cpld.etkf')
plt.plot(days,err3,label = 'true-solo.etkf')
plt.xlabel('days')
plt.ylabel('RMSE')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.title('ETKF; Ens = %s' % str(ens_num))
plt.savefig("rmse_etkf_ocn_%s_%3.1f%s.png" % (str(ens_num),infl,tw_da))

#plt.show()
