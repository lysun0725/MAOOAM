import numpy as np
import matplotlib.pyplot as plt
import sys

ens_num = sys.argv[1]
tw_da = sys.argv[2]
tw_solo = sys.argv[3]
tw_solo2 = sys.argv[4]
infl = 1.0
da_file = "Xam_etkf_cpl_%s_%3.1f%s_03.dat" % (str(ens_num),infl,tw_da)
wc_file = "Xam_etkf_wcp_%s_%3.1f%s_03.dat" % (str(ens_num),infl,tw_da)
sp_file = "../../MAOOAM_solo_atm/fortran_da/h%s/Xam_etkf_%s_%3.1f%s.dat" % (tw_solo,str(37),infl,tw_da)
sp_file2 = "../../MAOOAM_solo_ocn/fortran_da/h%s/Xam_etkf_%s_%3.1f%s.dat" % (tw_solo,str(37),infl,tw_da)
sp_file3 = "../../MAOOAM_solo_atm/fortran_da/h%s/Xam_etkf_%s_%3.1f%s.dat" % (tw_solo2,str(37),infl,tw_da)
sp_file4 = "../../MAOOAM_solo_ocn/fortran_da/h%s/Xam_etkf_%s_%3.1f%s.dat" % (tw_solo2,str(37),infl,tw_da)
fr_file = "../../MAOOAM_solo_atm/fortran_da/freerun_atm.010.dat"
fr_file2 = "../../MAOOAM_solo_ocn/fortran_da/freerun_ocn.010.dat"


true = np.loadtxt('nature.dat')[1:900002,1:37]
etkf = np.loadtxt(da_file)[:,1:37]
wcda = np.loadtxt(wc_file)[:,1:37]
sprd = np.loadtxt(sp_file)[:,1:21]
sprd2 = np.loadtxt(sp_file2)[0:900001,1:17]
sprd3 = np.loadtxt(sp_file3)[:,1:21]
sprd4 = np.loadtxt(sp_file4)[0:900001,1:17]
free = np.loadtxt(fr_file)[1:900002,1:21]
free2 = np.loadtxt(fr_file2)[1:900002,21:37]
days = np.arange(1,900002)*0.1/8.9

# atm
err6 = np.linalg.norm(true[:,0:20]-free,axis=1)/np.sqrt(20)
err2 = np.linalg.norm(true[:,0:20]-etkf[:,0:20],axis=1)/np.sqrt(20)
err8 = np.linalg.norm(true[:,0:20]-wcda[:,0:20],axis=1)/np.sqrt(20)
err3 = np.linalg.norm(true[:,0:20]-sprd,axis=1)/np.sqrt(20)
err10 = np.linalg.norm(true[:,0:20]-sprd3,axis=1)/np.sqrt(20)

# ocn
err7 = np.linalg.norm(true[:,20:36]-free2,axis=1)/np.sqrt(16)
err4 = np.linalg.norm(true[:,20:36]-etkf[:,20:36],axis=1)/np.sqrt(16)
err9 = np.linalg.norm(true[:,20:36]-wcda[:,20:36],axis=1)/np.sqrt(16)
err5 = np.linalg.norm(true[:,20:36]-sprd2,axis=1)/np.sqrt(16)
err11 = np.linalg.norm(true[:,20:36]-sprd4,axis=1)/np.sqrt(16)

plt.figure(1)
plt.plot(days,err6,label = 'freerun')
plt.plot(days,err2,label = 'cpld ETKF')
plt.plot(days,err8,label = 'wcpl ETKF')
plt.plot(days,err3,label = 'uncpld ETKF; h=%s' % tw_solo)
plt.plot(days,err10,label = 'uncpld ETKF; h=%s' % tw_solo2)
plt.xlabel('days')
plt.ylabel('RMSE')
plt.yscale('log')
plt.xscale('log')
plt.title('ETKF; Atm; Ens_num = %s' % str(ens_num))
plt.legend()

plt.savefig("rmse_etkf_catm_%s_%s_%3.1f%s_2.png" % (tw_solo,str(ens_num),infl,tw_da))

plt.figure(2)
plt.plot(days,err7,label = 'freerun')
plt.plot(days,err4,label = 'cpld ETKF')
plt.plot(days,err9,label = 'wcpl ETKF')
plt.plot(days,err5,label = 'uncpld ETKF; h=%s' % tw_solo)
plt.plot(days,err11,label = 'uncpld ETKF; h=%s' % tw_solo2)
plt.xlabel('days')
plt.ylabel('RMSE')
plt.yscale('log')
plt.xscale('log')
plt.title('ETKF; Ocn; Ens_num = %s' % str(ens_num))
plt.legend()

plt.savefig("rmse_etkf_cocn_%s_%s_%3.1f%s_2.png" % (tw_solo,str(ens_num),infl,tw_da))

#plt.show()
