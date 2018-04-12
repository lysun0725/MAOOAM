# coding: utf-8

import numpy as np
import tl_ad

def compute_tlm_solo(eta_atm,eta_ocn,dt,tensor,solo_flag):
    natm = np.shape(eta_atm)[0]
    noc = np.shape(eta_ocn)[0]
    ndim = natm + noc

    eta = np.concatenate((eta_atm,eta_ocn))
    tlm_cpld = tl_ad.compute_tlm(eta,dt,tensor)
    
    if solo_flag == 'atm':
        tlm_solo = tlm_cpld[0:natm,0:natm]
    elif solo_flag == 'ocn':
        tlm_solo = tlm_cpld[natm:ndim,natm:ndim]
    else:
        print ('ERROR...')

    return tlm_solo
    
def compute_adm_solo(eta_atm,eta_ocn,dt,tensor,solo_flag):
    natm = np.shape(eta_atm)[0]
    noc = np.shape(eta_ocn)[0]
    ndim = natm + noc

    eta = np.concatenate((eta_atm,eta_ocn))
    adm_cpld = tl_ad.compute_adm(eta,dt,tensor)
    
    if solo_flag == 'atm':
        adm_solo = adm_cpld[0:natm,0:natm]
    elif solo_flag == 'ocn':
        adm_solo = adm_cpld[natm:ndim,natm:ndim]
    else:
        print ('ERROR...')

    return tlm_solo
 
