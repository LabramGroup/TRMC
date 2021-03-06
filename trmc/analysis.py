import scipy.optimize
import pandas as pd
import numpy as np
import re



def offsettime(df, timebefore = 0, timeafter = None):
    """remove all data 'timebefore' before the max of the dataframe, then move the max to zero time"""
#     df_offset = pd.DataFrame(columns = df.columns)
    timemax = df[df.columns[0]].idxmax()
    time = df.index
    time1 = timemax-timebefore
    idx1 = time.get_loc(time1, method = 'nearest')

    if timeafter is None:
        idx2 = -1
    else:
        time2 = timemax+timeafter
        idx2 = time.get_loc(time2, method = 'nearest')

    df_cut = df.iloc[idx1:idx2]
    df_cut = df_cut.set_index(time[idx1:idx2] - timemax)
    return df_cut

def calc_K(f0,R0_norm,w, printparams = False, coupling = 'under'):   
    """Calculate the K value from the lorentzian fit constants""" 
    Q = f0/w
    t_rc = Q/(np.pi*f0)

    cav_w = 22.86e-3 #km(?)
    cav_h = 10.16e-3 #km
    cav_l = 76e-3 #mm  just came up with this number

    beta = cav_w/cav_h
    eps_r = 1 #dielectric
    eps_0 = 8.85e-12
    
    if coupling == 'under':
        K = -( 2*Q*(( (1/np.sqrt(R0_norm)) - 1 ) ) )/(np.pi*f0*eps_r*eps_0*cav_l*beta)
    elif coupling == 'over':
        K = ( 2*Q*(( (1/np.sqrt(R0_norm)) + 1 ) ) )/(np.pi*f0*eps_r*eps_0*cav_l*beta)
    if(printparams):
#         print('eps: ', "{:.2E}".format(eps))
#         print('beta: ', "{:.2E}".format(beta))
#         print('cav_l: ', "{:.2E}".format(cav_l))
        print('f0: ', "{:.2E}".format(f0))
        print('w: ', "{:.2E}".format(w))
        print('R0: ', "{:.2E}".format(R0))
        print('Checks:')
        print('Q: ', "{:.2E}".format(Q))
        print('t_rc: ', "{:.2E}".format(t_rc))
        print('Output: ')
        print('K: ', "{:.2E}".format(K))
    return K

import itertools
def dict_product(dicts):
    """
    >>> list(dict_product(dict(number=[1,2], character='ab')))
    [{'character': 'a', 'number': 1},
     {'character': 'a', 'number': 2},
     {'character': 'b', 'number': 1},
     {'character': 'b', 'number': 2}]
    """
    return (dict(zip(dicts, x)) for x in itertools.product(*dicts.values()))


###Fitting###

def fitsweep(v, p0, bounds, window, fittype, p_labels):
    xdata = v.indexes['freq'].values
    ydata = v.values

    if fittype == 'lor':
        if p0[0] == None:
            f0 = v.to_series().idxmin()
            p0[0] = f0

        v_p, v_sl = fit_lor(xdata,ydata, p0,bounds, window)
        popt = v_p[0]
        v_fit = lorfn(*popt)

    if fittype == 'lorline':
        if p0[0] == None:
            f0 = v.to_series().idxmin()
            p0[0] = f0

        v_p, v_sl = fit_lor_line(xdata,ydata, p0,bounds, window)
        popt = v_p[0]
        v_fit = lor_line_fn(*popt)

    elif fittype == 'poly':
        
        v_p, v_sl = fit_poly2(xdata,ydata, p0, bounds , window)
        popt = v_p[0]
        popc = v_p[1]
        v_fit = poly2fn(*popt)
        minf, minR  = polymin(v, v_sl, v_fit)
        v_p = [[minf,minR, *popt],popc]

    return v_fit, v_p, v_sl


def fit_lor_line(xdata,ydata, p0, bounds = ([0,0,0, 0, -np.inf,0],[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]), window = 105):
    """Fits to lorentzian function and returns parameters"""
    #xdata = sweep.index.values
    #ydata = sweep.values

    minidx = ydata.argmin()
    minfreq = xdata[minidx]

    sl = slice(minidx-window,minidx+window+1)
    popt,popc = scipy.optimize.curve_fit(lor_line,xdata[sl],ydata[sl], p0 , bounds = bounds)
    p = (popt,popc)
    return p, sl

def lor_line(f,f0,R0,Rinf, w, m, b): 
    return (R0 + Rinf*(2*(f-f0)/w)**2)/(1 + (2*(f-f0)/w)**2) + m*(f - f0) + b

def lor_line_fn(f0,R0,Rinf,w ,m,b): 
    """retruns a lorentzian function defined by parameters"""
    def fn(f):
        return lor_line(f,f0,R0, Rinf,w,m,b)
    return fn 

def fit_lor(xdata,ydata, p0, bounds = ([0,0,0, 0],[np.inf,np.inf,np.inf,np.inf]), window = 105):
    """Fits to lorentzian function and returns parameters"""
    #xdata = sweep.index.values
    #ydata = sweep.values

    minidx = ydata.argmin()
    minfreq = xdata[minidx]

    sl = slice(minidx-window,minidx+window+1)
    popt,popc = scipy.optimize.curve_fit(lor,xdata[sl],ydata[sl], p0 , bounds = bounds)
    p = (popt,popc)
    return p, sl

def lorfn(f0,R0,Rinf,w ): 
    """retruns a lorentzian function defined by parameters"""
    def fn(f):
        return lor(f,f0,R0, Rinf,w)
    return fn 

def lor(f,f0,R0,Rinf, w): 
    return (R0 + Rinf*(2*(f-f0)/w)**2)/(1 + (2*(f-f0)/w)**2)


def fit_poly2(xdata,ydata, p0, bounds = ([0,0,0],[np.inf,np.inf,np.inf]), window = 105):
    """Fits to a polynomial and returns fit function and parameters"""
    # xdata = sweep.index.values
    # ydata = sweep.values

    minidx = ydata.argmin()
    # minfreq = xdata[minidx]

    sl = slice(minidx-window,minidx+window)

    popt,popc = scipy.optimize.curve_fit(poly2,xdata[sl],ydata[sl], p0 , bounds = bounds)
    p = (popt,popc)
    return p, sl

def poly2fn(c0,c1,c2):
    def fn(x):
        return poly2(x,c0,c1,c2)
    return fn

def poly2(x,c0,c1,c2):
    return c0 + c1*x + c2*x**2
    


def fit_poly_np(sweep, window = 105, order = 3):
    """Fits to a polynomial and returns fit function and parameters"""
    xdata = sweep.index.values
    ydata = sweep.values

    minidx = ydata.argmin()
    minfreq = xdata[minidx]

    bounds = ([0,0,0],[np.inf,np.inf,np.inf])

    sl = slice(minidx-window,minidx+window)
    p = np.polyfit(xdata[sl],ydata[sl], order)
    fit_func = np.poly1d(p)

    return fit_func, p, sl


def polymin(v0,v0_sl, fit):
    """Returns the minimum of a fit functon within the window defined by v0[v0_sl]"""
    f = np.linspace(v0.indexes['freq'][v0_sl][0],v0.indexes['freq'][v0_sl][-1],num = 1000)
    fitdata = fit(f)
    minR = fitdata.min()
    minf = f[fitdata.argmin()]
    return minf,minR




###CONVENTIONAL ANALYSIS###


def convert_V2cond(df_V,back_V,K):
    df_cond = - ((df_V)/back_V)/K
    return df_cond

def maxG_and_fom(df_cond, params):
    """Calculates maxG and the figure of merit from a dataframe of the deltaG"""
    beta = params['beta']
    e = 1.6e-19
    FA = params['FA']*0.9 #0.9 factor from ITO
    M = params['M']

    fluences = df_cond.columns
    maxG = pd.Series(index = fluences)
    fom = pd.Series(index = fluences)

    for fluence in fluences:
        maxG[fluence] = df_cond[fluence].max()
        fom[fluence] =  maxG[fluence]/(beta*e*fluence*FA*M)

    return maxG, fom
