import pandas as pd
import numpy as np


def calc_pow(t,I0,params):
    """Calulcates the time dependent power of the pulse"""
    sig = params['FWHM']/(2*np.sqrt(2*np.log(2)))
    p0 = (I0/(sig*np.sqrt(2*np.pi)))
    power = p0*np.exp(-((t-params['t0'])**2)/(2*(sig**2)))
    power = pd.Series(power, index = t)
    dng = (params['FA']/params['d'])*power
    return power, dng

def dnr(n,k1,k2,k3):
    """Recombination rate"""
    dnr = np.zeros(len(n))
    if(k1 != 0):
        dnr = dnr - k1*n
    if(k2 != 0):
        dnr = dnr - k2*n**2
    if(k3 != 0):
        dnr = dnr - k3*n**3
    return dnr

def calc_n(dng,k1,k2,k3):
    """numerical integration to find number density"""
    t = dng.index
    n = pd.Series(np.zeros(len(t)),index = t)
    for i in range(1,len(t)):
        dng_t = dng.iloc[:i-1]
        dnr_t = dnr(n.iloc[:i-1],k1,k2,k3)
        dn_t = dng_t + dnr_t
        n.iloc[i] = np.trapz(dn_t,t[:i-1])
    return n