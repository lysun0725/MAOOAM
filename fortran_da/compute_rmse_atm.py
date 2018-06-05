import numpy as np
import matplotlib.pyplot as plt
import sys

ens_num = sys.argv[1]
tw_da = sys.argv[2]
tw_da2 = ''
infl = 1.0
da_file = "Xam_etkf_atm_%s_%3.1f%s.dat" % (str(ens_num),infl,tw_da)
#sp_file = "sprd_etkf_ocn_%s_%3.1f%s.dat" % (str(ens_num),infl,tw_da)
sp_file = "../../MAOOAM_solo_atm/fortran_da/h.1E-01/Xam_etkf_%s_%3.1f%s.dat" % (str(20),infl,tw_da2)

true = np.loadtxt('nature.dat')[1:900002,1:21]
#free = np.loadtxt('freerun_atm.010.dat')[1:900002,1:21]
etkf = np.loadtxt(da_file)[:,1:21]
sprd = np.loadtxt(sp_file)[:,1:21]

#err1 = np.linalg.norm(true-free,axis=1)/np.sqrt(36)
err2 = np.linalg.norm(true-etkf,axis=1)/np.sqrt(20)
err3 = np.linalg.norm(true-sprd,axis=1)/np.sqrt(20)

#plt.plot(err1,label = 'true-solo')
plt.plot(err2,label = 'true-cpld.etkf')
plt.plot(err3,label = 'true-solo.etkf')
plt.xlabel('DA cycle')
plt.ylabel('RMSE')
plt.yscale('log')
plt.xscale('log')
plt.title('ETKF; Ens = %s' % str(ens_num))
plt.legend()

plt.savefig("rmse_etkf_atm_%s_%3.1f%s.png" % (str(ens_num),infl,tw_da))

#plt.show()
