# This code is for a plot of drifter trajectory vs. ocean stream function field.

import numpy as np
import matplotlib.pyplot as plt

def base_o(xp,yp,Ho,Po,n):
    res = 2*np.sin(Ho*n/2*xp)*np.sin(Po*yp) # eqn (9) in De Cruz et al. 2016
    return res

def base_a_A(yp,P):
    res = np.sqrt(2)*np.cos(P*yp) # eqn (6) in De Cruz et al. 2016
    return res

def base_a_K(xp,yp,M,P,n):
    res = 2*np.cos(M*n*xp)*np.sin(P*xp) # eqn (7) in De Cruz et al. 2016
    return res

def base_a_L(xp,yp,H,P,n):
    res = 2*np.sin(H*n*xp)*np.sin(P*yp) # eqn (8) in De Cruz et al. 2016
    return res

infile = 'd001/nature.dat'

atm_str = np.loadtxt(infile)[:,1:11] # atmosphere stream function nodes
atm_the = np.loadtxt(infile)[:,11:21] # atmosphere temp function nodes
ocn_str = np.loadtxt(infile)[:,21:29] # ocean stream function nodes
ocn_the = np.loadtxt(infile)[:,29:37] # ocean temp function nodes
dr_xy   = np.loadtxt(infile)[:,37:]

n = 1.5 # in params.nml
Hmax = 2; Pmax = 2 # De Cruz et al. 2016, Page 3
Ho = 2; Po = 4 # in params
Na = Pmax*(2*Hmax+1); No = Po*Ho
Nxp = 400; Nyp = 300 # depends on n
Nt = atm_str.shape[0]
dr_num = 1; dr_size = 2
dr_xy = np.reshape(dr_xy,[Nt,dr_num,dr_size])
xpvec = np.linspace(0,2*np.pi/n,Nxp)
ypvec = np.linspace(0,np.pi,Nyp)
xp, yp = np.meshgrid(xpvec,ypvec)
psi_a = np.zeros([Nyp,Nxp])
the_a = np.zeros([Nyp,Nxp])
psi_o = np.zeros([Nyp,Nxp]) # Define the terminal stream function field
the_o = np.zeros([Nyp,Nxp]) # Define the terminal stream function field
# the_o = np.zeros([Nxp,Nyp])

# compute atmos stream function field: info from IC.nml and De Cruz et al. 2016
# IC(1): typ=A; P=1    
# IC(2): typ=K; M=1 P=1
# IC(3): typ=L; H=1 P=1
# IC(4): typ=A; P=2    
# IC(5): typ=K; M=2 P=1
# IC(6): typ=L; H=2 P=1
# IC(7): typ=K; M=1 P=2
# IC(8): typ=L; H=1 P=2
# IC(9): typ=K; M=2 P=2
# IC(10):typ=L; H=2 P=2
psi_a = atm_str[-1,0]*base_a_A(yp,1) + atm_str[-1,3]*base_a_A(yp,2) 
psi_a = psi_a + atm_str[-1,1]*base_a_K(xp,yp,1,1,n) + atm_str[-1,2]*base_a_L(xp,yp,1,1,n)
psi_a = psi_a + atm_str[-1,4]*base_a_K(xp,yp,2,1,n) + atm_str[-1,5]*base_a_L(xp,yp,2,1,n)
psi_a = psi_a + atm_str[-1,6]*base_a_K(xp,yp,1,2,n) + atm_str[-1,7]*base_a_L(xp,yp,1,2,n)
psi_a = psi_a + atm_str[-1,8]*base_a_K(xp,yp,2,2,n) + atm_str[-1,9]*base_a_L(xp,yp,2,2,n)

# compute atmos temperature
the_a = atm_the[-1,0]*base_a_A(yp,1) + atm_the[-1,3]*base_a_A(yp,2)
the_a = the_a + atm_the[-1,1]*base_a_K(xp,yp,1,1,n) + atm_the[-1,2]*base_a_L(xp,yp,1,1,n)
the_a = the_a + atm_the[-1,4]*base_a_K(xp,yp,2,1,n) + atm_the[-1,5]*base_a_L(xp,yp,2,1,n)
the_a = the_a + atm_the[-1,6]*base_a_K(xp,yp,1,2,n) + atm_the[-1,7]*base_a_L(xp,yp,1,2,n)
the_a = the_a + atm_the[-1,8]*base_a_K(xp,yp,2,2,n) + atm_the[-1,9]*base_a_L(xp,yp,2,2,n)


# compute ocean stream function field.
cnt = 0
for i in np.arange(1,Ho+1):
  for j in np.arange(1,Po+1):
    const = 2*((-1)**i-1.0)*((-1)**j-1.0)/i/j/np.pi/np.pi # eqn (15) in De Cruz et al.2016
    psi_o = psi_o + ocn_str[-1,cnt] * (base_o(xp,yp,i,j,n) - const) # eqn (13) in  De Cruz et al. 2016 
    cnt += 1

# compute ocean temperature field.
cnt = 0
for i in np.arange(1,Ho+1):
  for j in np.arange(1,Po+1):
    the_o = the_o + ocn_the[-1,cnt] * base_o(xp,yp,i,j,n)  # eqn (14) in  De Cruz et al. 2016 
    cnt += 1

figure, axes = plt.subplots(nrows=2, ncols=2)
plt.subplot(2,2,1)
plt.contourf(xp,yp,psi_a,cmap='rainbow')
plt.colorbar()
plt.title('atmos psi')

plt.subplot(2,2,2)
plt.contourf(xp,yp,the_a,cmap='rainbow')
plt.colorbar()
plt.title('atmos T')

plt.subplot(2,2,3)
ax=plt.contourf(xp,yp,psi_o,cmap='rainbow')
plt.colorbar(ax,format='%.0e')
plt.plot(dr_xy[:,0,0],dr_xy[:,0,1])
plt.title('ocean psi')

plt.subplot(2,2,4)
plt.contourf(xp,yp,the_o,cmap='rainbow')
plt.colorbar()
plt.title('ocean theta')
figure.tight_layout()

plt.savefig('test_dr.eps')
plt.show()
