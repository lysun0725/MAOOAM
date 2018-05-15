import numpy as np
import matplotlib.pyplot as plt
import sys

ens_num = sys.argv[1]
tw_da = sys.argv[2]
infl = 1.0
da_file = "Xam_etkf_%s_%3.1f%s.dat" % (str(ens_num),infl,tw_da)
sp_file = "sprd_etkf_%s_%3.1f%s.dat" % (str(ens_num),infl,tw_da)

true = np.loadtxt('nature.dat')[1:900002,1:21]
free = np.loadtxt('freerun_atm.010.dat')[1:900002,1:21]
etkf = np.loadtxt(da_file)[:,1:21]
sprd = np.loadtxt(sp_file)[:,1:21]

err1 = np.linalg.norm(true-free,axis=1)/np.sqrt(20)
err2 = np.linalg.norm(true-etkf,axis=1)/np.sqrt(20)
err3 = np.linalg.norm(sprd,axis=1)/np.sqrt(20)

plt.plot(err1,label = 'true-free')
plt.plot(err2,label = 'true-etkf')
plt.plot(err3,label = 'sprd')
plt.xlabel('DA cycle')
plt.ylabel('RMSE')
plt.yscale('log')
plt.xscale('log')
plt.legend()

plt.savefig("rmse_etkf_%s_%3.1f%s.png" % (str(ens_num),infl,tw_da))

#plt.show()
