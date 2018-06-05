import numpy as np
import matplotlib.pyplot as plt
import sys

ens_num = sys.argv[1]
tw_da = sys.argv[2]
ens_num2 = sys.argv[3]
infl = 1.0
da_file = "Xam_etkf_cpl_%s_%3.1f%s.dat" % (str(ens_num),infl,tw_da)
sp_file = "Xam_etkf_atm_%s_%3.1f%s.dat" % (str(ens_num2),infl,tw_da)
sp_file2 = "Xam_etkf_ocn_%s_%3.1f%s.dat" % (str(ens_num2),infl,tw_da)
sp_file3 = "sprd_etkf_cpl_%s_%3.1f%s.dat" % (str(ens_num),infl,tw_da)
sp_file4 = "sprd_etkf_atm_%s_%3.1f%s.dat" % (str(ens_num),infl,tw_da)
sp_file5 = "sprd_etkf_ocn_%s_%3.1f%s.dat" % (str(ens_num),infl,tw_da)
fr_file = "freerun_atm.010.dat"
fr_file2 = "freerun_ocn.010.dat"

true = np.loadtxt('nature.dat')[1:900002,1:37]
etkf = np.loadtxt(da_file)[:,1:37]
sprd = np.loadtxt(sp_file)[:,1:37]
sprd2= np.loadtxt(sp_file2)[:,1:37]
sprd3= np.loadtxt(sp_file3)[:,1:37]
sprd4= np.loadtxt(sp_file4)[:,1:37]
sprd5= np.loadtxt(sp_file5)[:,1:37]
free = np.loadtxt(fr_file)[1:900002,1:21]
free2 = np.loadtxt(fr_file2)[1:900002,21:37]
days = np.arange(1,900002)*0.1/8.9

# atm
err2 = np.linalg.norm(true[:,0:20]-etkf[:,0:20],axis=1)/np.sqrt(20)
err3 = np.linalg.norm(true[:,0:20]-sprd[:,0:20],axis=1)/np.sqrt(20)
err1 = np.linalg.norm(true[:,0:20]-sprd2[:,0:20],axis=1)/np.sqrt(20)
err7 = np.linalg.norm(sprd3[:,0:20],axis=1)/np.sqrt(20)
err8 = np.linalg.norm(sprd4[:,0:20],axis=1)/np.sqrt(20)
err9 = np.linalg.norm(sprd5[:,0:20],axis=1)/np.sqrt(20)
err13 = np.linalg.norm(true[:,0:20]-free,axis=1)/np.sqrt(20)

# ocn
err4 = np.linalg.norm(true[:,20:36]-etkf[:,20:36],axis=1)/np.sqrt(16)
err5 = np.linalg.norm(true[:,20:36]-sprd[:,20:36],axis=1)/np.sqrt(16)
err6 = np.linalg.norm(true[:,20:36]-sprd2[:,20:36],axis=1)/np.sqrt(16)
err10 = np.linalg.norm(sprd3[:,20:36],axis=1)/np.sqrt(16)
err11 = np.linalg.norm(sprd4[:,20:36],axis=1)/np.sqrt(16)
err12 = np.linalg.norm(sprd5[:,20:36],axis=1)/np.sqrt(16)
err14 = np.linalg.norm(true[:,20:36]-free2,axis=1)/np.sqrt(16)

plt.figure(1)
plt.plot(days,err13,label = 'freerun')
plt.plot(days,err2,label = 'obs:cpld; ens:%s' % str(ens_num))
plt.plot(days,err3,label = 'obs:atm; ens:%s' % str(ens_num2))
plt.plot(days,err1,label = 'obs:ocn; ens:%s' % str(ens_num2))
plt.xlabel('days')
plt.ylabel('RMSE')
plt.yscale('log')
plt.xscale('log')
plt.title('RMSE for atmosphere components; tw_run = %s' % str(tw_da))
plt.legend()

plt.savefig("rmse_etkf_catm_%s_%s_%3.1f%s.png" % (str(ens_num),str(ens_num2),infl,tw_da))

plt.figure(3)
plt.plot(days,err13,label = 'freerun')
plt.plot(days,err7,label = 'obs:cpld; ens:%s' % str(ens_num))
plt.plot(days,err8,label = 'obs:atm; ens:%s' % str(ens_num2))
plt.plot(days,err9,label = 'obs:ocn; ens:%s' % str(ens_num2))
plt.xlabel('days')
plt.ylabel('Ensemble spreads')
plt.yscale('log')
plt.xscale('log')
plt.title('Ensemble spreads for atmosphere components; tw_run = %s' % str(tw_da))
plt.legend()

plt.savefig("rmse_etkf_catm_%s_%s_%3.1f%s_sprd.png" % (str(ens_num),str(ens_num2),infl,tw_da))

plt.figure(2)
plt.plot(days,err14,label = 'freerun')
plt.plot(days,err4,label = 'obs:cpld; ens:%s' % str(ens_num))
plt.plot(days,err5,label = 'obs:atm; ens:%s' % str(ens_num2))
plt.plot(days,err6,label = 'obs:ocn; ens:%s' % str(ens_num2))
plt.xlabel('days')
plt.ylabel('RMSE')
plt.yscale('log')
plt.xscale('log')
plt.title('RMSE for ocean components; tw_run = %s' % str(tw_da))
plt.legend()

plt.savefig("rmse_etkf_cocn_%s_%s_%3.1f%s.png" % (str(ens_num),str(ens_num2),infl,tw_da))

plt.figure(4)
plt.plot(days,err14,label = 'freerun')
plt.plot(days,err10,label = 'obs:cpld; ens:%s' % str(ens_num))
plt.plot(days,err11,label = 'obs:atm; ens:%s' % str(ens_num2))
plt.plot(days,err12,label = 'obs:ocn; ens:%s' % str(ens_num2))
plt.xlabel('days')
plt.ylabel('Ensemble spreads')
plt.yscale('log')
plt.xscale('log')
plt.title('Ensemble spreads for ocean components; tw_run = %s' % str(tw_da))
plt.legend()

plt.savefig("rmse_etkf_cocn_%s_%s_%3.1f%s_sprd.png" % (str(ens_num),str(ens_num2),infl,tw_da))

#plt.show()
