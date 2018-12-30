import pandas as pd
import numpy as np
import re


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

def lorfn(f0,w,R0, Rinf): 
    """retruns a lorentzian function defined by parameters"""
    def fn(f):
        return lor(f,f0,w,R0, Rinf)
    return fn 

def lor(f,f0,w,R0, Rinf): 
    return (R0 + Rinf*(2*(f-f0)/w)**2)/(1 + (2*(f-f0)/w)**2)


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


def calc_K(f0,w,R0, printparams = False):   
    """Calculate the K value from the lorentzian fit constants""" 
    Q = f0/w
    t_rc = Q/(np.pi*f0)

    cav_w = 22.86e-3 #km(?)
    cav_h = 10.16e-3 #km
    cav_l = 76e-3 #mm  just came up with this number

    beta = cav_w/cav_h
    eps_r = 1 #dielectric
    eps_0 = 8.85e-12
    
    K = ( 2*Q*( 1 + (1/np.sqrt(R0)) ) )/(np.pi*f0*eps_r*eps_0*cav_l*beta)
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

if __name__ == '__main__':
    filepaths_A = ['C:\\Users\\aspit\\OneDrive\\Data\\TRMC\\Gratzel\\Sample A\\Data\\High_Power_Filter=01_Fluence=6.45E+14_data.csv','C:\\Users\\aspit\\OneDrive\\Data\\TRMC\\Gratzel\\Sample A\\Data\\High_Power_Filter=02_Fluence=5.121E+14_data.csv', 'C:\\Users\\aspit\\OneDrive\\Data\\TRMC\\Gratzel\\Sample A\\Data\\High_Power_Filter=03_Fluence=4.07E+14_data.csv', 'C:\\Users\\aspit\\OneDrive\\Data\\TRMC\\Gratzel\\Sample A\\Data\\High_Power_Filter=04_Fluence=3.231E+14_data.csv', 'C:\\Users\\aspit\\OneDrive\\Data\\TRMC\\Gratzel\\Sample A\\Data\\High_Power_Filter=05_Fluence=2.567E+14_data.csv', 'C:\\Users\\aspit\\OneDrive\\Data\\TRMC\\Gratzel\\Sample A\\Data\\High_Power_Filter=06_Fluence=2.038E+14_data.csv', 'C:\\Users\\aspit\\OneDrive\\Data\\TRMC\\Gratzel\\Sample A\\Data\\High_Power_Filter=07_Fluence=1.619E+14_data.csv', 'C:\\Users\\aspit\\OneDrive\\Data\\TRMC\\Gratzel\\Sample A\\Data\\High_Power_Filter=08_Fluence=6.45E+13_data.csv', 'C:\\Users\\aspit\\OneDrive\\Data\\TRMC\\Gratzel\\Sample A\\Data\\High_Power_Filter=09_Fluence=3.231E+13_data.csv', 'C:\\Users\\aspit\\OneDrive\\Data\\TRMC\\Gratzel\\Sample A\\Data\\High_Power_Filter=10_Fluence=6.45E+12_data.csv', 'C:\\Users\\aspit\\OneDrive\\Data\\TRMC\\Gratzel\\Sample A\\Data\\High_Power_Filter=11_Fluence=6.45E+11_data.csv']
    df_A_V, df_A_cond = load(filepaths_A, offsettime = 50e-9, sub_lowpow = True)
    print(df_A_V)
