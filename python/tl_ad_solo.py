# coding: utf-8

import numpy as np
import tl_ad
from integrator import tendencies

def jacobi_solo_mat(eta_atm,eta_ocn,tensor,solo_flag):
    natm = np.shape(eta_atm)[0]
    noc = np.shape(eta_ocn)[0]
    ndim = natm + noc
    
    if solo_flag == 'atm':
        j_mat = np.zeros([natm,natm])
    elif solo_flag == 'ocn':
        j_mat = np.zeros([noc,noc])

    eta2 = np.append(np.append(eta_atm, eta_ocn),1)
        
    for n in np.arange(0,np.shape(tensor)[0]):
        i = tensor[n][0]
        j = tensor[n][1]
        k = tensor[n][2]
        v = tensor[n][3]

        if solo_flag == 'atm':  # for solo_atm
            if j!=0 and i<=natm and j<=natm:
                j_mat[i-1,j-1] = j_mat[i-1,j-1] + v*eta2[k-1]
        
            if k!=0 and i<=natm and k<=natm:
                j_mat[i-1,k-1] = j_mat[i-1,k-1] + v*eta2[j-1]
        elif solo_flag == 'ocn': # for solo_ocn
            if j!=0 and i>natm and j>natm:
                j_mat[i-1-natm,j-1-natm] = j_mat[i-1-natm,j-1-natm] + v*eta2[k-1]
        
            if k!=0 and i>natm and k>natm:
                j_mat[i-1-natm,k-1-natm] = j_mat[i-1-natm,k-1-natm] + v*eta2[j-1]
        else:
            print ('ERROR...')

    return j_mat

def compute_tlm_solo(eta_atm,eta_ocn,dt,tensor,solo_flag):
    natm = np.shape(eta_atm)[0]
    noc = np.shape(eta_ocn)[0]
    ndim = natm + noc

    #j_mat = jacobi_solo_mat(eta_atm,eta_ocn,tensor,solo_flag)
    if solo_flag == 'atm':
        j1 = jacobi_solo_mat(eta_atm,eta_ocn,tensor,solo_flag)
        eta = np.append(eta_atm,eta_ocn)
        k1 = tendencies(eta)[0:natm]
        eta_atm1 = eta_atm + k1*dt
        j2 = jacobi_solo_mat(eta_atm1,eta_ocn,tensor,solo_flag)
        tlm_solo = np.identity(natm) + (j1+j2)*0.5*dt + np.dot(j2,j1)*0.5*dt**2
    elif solo_flag == 'ocn':
        j1 = jacobi_solo_mat(eta_atm,eta_ocn,tensor,solo_flag)
        eta = np.append(eta_atm,eta_ocn)
        k1 = tendencies(eta)[natm:ndim]
        eta_ocn1 = eta_ocn + k1*dt
        j2 = jacobi_solo_mat(eta_atm,eta_ocn1,tensor,solo_flag)
        tlm_solo = np.identity(noc) + (j1+j2)*0.5*dt + np.dot(j2,j1)*0.5*dt**2

    else:
        print ('ERROR...')

    return tlm_solo
    
def compute_adm_solo(eta_atm,eta_ocn,dt,tensor,solo_flag):
    natm = np.shape(eta_atm)[0]
    noc = np.shape(eta_ocn)[0]
    ndim = natm + noc

    j_mat = jacobi_solo_mat(eta_atm,eta_ocn,tensor,solo_flag)
    if solo_flag == 'atm':
        adm_solo = np.identity(natm) + dt*j_mat.transpose()
    elif solo_flag == 'ocn':
        adm_solo = np.identity(noc) + dt*j_mat.transpose()
    else:
        print ('ERROR...')

    return adm_solo
 
