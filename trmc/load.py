import pandas as pd
import numpy as np
import re


def read_params(filepath):
    """Read the amplification and background voltage from the header of a file"""
    params = pd.read_csv(filepath, nrows = 11, usecols = [1])
    params = params.transpose()
    amp = float(params['Amplification'][0]) # Is amplification already taken into account?
    # K = float(params['K'][0])
    back_V = params['Background Voltage'][0].replace('V','')
    unitdict = {'m':1e-3, 'u':1e-6}
    scale = unitdict[back_V[-1]]
    back_V = float(back_V[:len(back_V)-1])*scale # remove mV or uV and scale appropriately
    return amp, back_V

def load_trace(filepath,offsettime = None):
    """load in a single trace csv file"""
    temp = pd.read_csv(filepath, skiprows = 13, index_col=0)
    volt = temp['Voltage (V)']
    if offsettime is not None:
        volt = volt - np.mean(volt[0:offsettime])

    return volt
    # time = volt.index

def load_fluencesweep(filepaths, offsettime = None, sub_lowpow = False):
    """
    Load in a csv set of flucence sweep data. This sets columns to fluence values which is not tidy data and should be changed.

    offsettime - takes an average of the data between 0 and offsettime and subtracts that from the data
    sub_lowpow - subtracts the low power trace from all datasets
    """

    V1 = pd.read_csv(filepaths[-1], skiprows = 13,index_col = 0)['Voltage (V)']
    time  = V1.index
    df_V= pd.DataFrame(index = time)
    fluences = pd.Series(index = filepaths)
    for filepath in filepaths:
        volt = read_trace(filepath, offsettime=offsettime)
        if(sub_lowpow):
            volt = volt - V1
        df_V = pd.concat([df_V, volt], axis = 1)
        
        m = re.search('Fluence=(.+?)_',filepath) # Read fluence from filename
        if m:
            fluences[filepath] = m.group(1)
    df_V.columns = fluences

    if(sub_lowpow):
        df_V = df_V.drop(df_V.columns[-1], axis = 1)

    # df_cond = convert_V2cond(df_V,back_V,K)

    return df_V



def load_fluence(filepath):
    """Loads in fluence data, not sure if working"""
    # filepath = 'C:\\Users\\aspit\\OneDrive\\Data\\TRMC\\Gratzel\\Sample A\\High_Power_3_FluenceSweep.csv'
    fluencesweep = pd.read_csv(filepath, skiprows = 13)
    fluences = fluencesweep['Fluence(cm^-2)']
    return fluences
