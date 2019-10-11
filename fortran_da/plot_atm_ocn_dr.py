import numpy as np
import matplotlib.pyplot as plt

def base_o(xp,yp,Ho,Po,n):
    res = 2*np.sin(Ho*n/2*xp)*np.sin(Po*yp)
    return res

atm_str = np.loadtxt('evol_field_dr.dat')[:,1:11] # atmosphere stream function nodes
atm_the = np.loadtxt('evol_field_dr.dat')[:,11:21] # atmosphere temp function nodes
ocn_str = np.loadtxt('evol_field_dr.dat')[:,21:29] # ocean stream function nodes
ocn_the = np.loadtxt('evol_field_dr.dat')[:,29:37] # ocean temp function nodes
dr_xy   = np.loadtxt('evol_field_dr.dat')[:,37:]

n = 1.5 # in params.nml
Ho = 2; Po = 4 # in params
Nxp = 400; Nyp = 300 # depends on n
Nt = atm_str.shape[0]
dr_num = 1; dr_size = 2
dr_xy = np.reshape(dr_xy,[Nt,dr_num,dr_size])
xpvec = np.linspace(0,2*np.pi/n,Nxp)
ypvec = np.linspace(0,np.pi,Nyp)
xp, yp = np.meshgrid(xpvec,ypvec)
psi_o = np.zeros([Nyp,Nxp])
# the_o = np.zeros([Nxp,Nyp])

# compute ocean stream function field.
cnt = 0
for i in np.arange(1,Ho+1):
  for j in np.arange(1,Po+1):
    const = 2*((-1)**i-1.0)*((-1)**j-1.0)/i/j/np.pi/np.pi
    psi_o = psi_o + ocn_str[-1,cnt] * (base_o(xp,yp,i,j,n) - const)      
    cnt += 1

plt.contourf(xp,yp,psi_o,cmap='rainbow')
plt.colorbar()
plt.plot(dr_xy[:,0,0],dr_xy[:,0,1])
plt.savefig('test.eps')
