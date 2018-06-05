
# coding: utf-8

# In[24]:


import numpy as np
import aotensor as aotensor_mod
import integrator

def jacobi_mat(eta,tensor):
    ndim = np.shape(eta)[0]
    j_mat = np.zeros([ndim,ndim])
   
    eta2 = np.append(eta,1.) 
    for n in np.arange(0,np.shape(tensor)[0]):
        i = tensor[n][0]
        j = tensor[n][1]
        k = tensor[n][2]
        v = tensor[n][3]
    
        if j!=0:
            j_mat[i-1,j-1] = j_mat[i-1,j-1] + v*eta2[k-1]
        
        if k!=0:
            j_mat[i-1,k-1] = j_mat[i-1,k-1] + v*eta2[j-1]

    return j_mat

def compute_tlm(eta,dt,tensor):
    ndim = np.shape(eta)[0]
    n = 1
    eta1 = eta
    j1 = jacobi_mat(eta1,tensor)
    eta2 = integrator.step(eta1,0,dt)
    j2 = jacobi_mat(eta2,tensor)
    
    tlm = np.identity(ndim) + (j1+j2)*0.5*dt + np.dot(j2,j1)*0.5*dt**2
    
    return tlm

def tlm_validation(eta,dt,tlm):
    ndim = np.shape(eta)[0]
    err = 0.1*np.random.randn(ndim)
    eta_new = integrator.step(eta,0,dt)
    eta_new_p = integrator.step(eta+err,0,dt)

    out = np.linalg.norm(eta_new-eta_new_p)/np.linalg.norm(tlm*err)

    return out

def compute_adm(eta,dt,tensor):

    # Notice, the input dt is a positive number
    ndim = np.shape(eta)[0]
    n = 1
    dt_adm = dt/n
    eta_old =  eta
    adm_old = np.identity(ndim)

    for i in range(1,n+1):
        j_mat = jacobi_mat(eta_old,tensor)
        adm = np.dot(np.identity(ndim)+dt_adm*j_mat.transpose(),adm_old)
       
        eta_new = integrator.step(eta_old,0,-dt_adm)
        eta_old = eta_new 
        adm_old = adm
    return adm


# In[26]:





    
    
        

