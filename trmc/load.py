import os
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


def freqfluence_flist(direc):
    """Creates a multindexed Series of filepaths from a frequency fluence sweep folder"""
    folders = os.listdir(direc)
    miarray = []
    folder_re = '^(\d+\.\d+)GHz_(.+?)'
    file_re = 'FF_Filter=\d+_Fluence=(.+?)_data.csv'

    flist = []

    for folder in folders:
        m_folder = re.search(folder_re,folder)
        freq = float(m_folder.groups(0)[0])*1e9
        direction = m_folder.groups(0)[1]
        folderpath = os.path.join(direc,folder)
        files = os.listdir(folderpath)
        for file in files:
            if file[0] == 'F':
                m_file = re.search(file_re, file)
                fluence = float(m_file.groups(0)[0])

                fp = os.path.join(folderpath,file)

                miarray.append((direction,freq,fluence))
                flist.append(fp)


    mi = pd.MultiIndex.from_tuples(miarray)

    s_fps = pd.Series(flist, index = mi)
    
    return s_fps


def freqfluence_load(s_fps):
    """Takes in a frequecny fluence sweep filepath Series and loads into a data Series"""
    direcs = set(s_fps.index.levels[0])
    freqs = set(s_fps.index.levels[1])
    fluences = set(s_fps.index.levels[2])
    
    mi = s_fps.index
    time_arr = load_trace(s_fps.iloc[0]).index.values
    miarray_t = []
    for tup in mi:
        for time in time_arr:
            miarray_t.append((*tup,time))

    mi_t = pd.MultiIndex.from_tuples(miarray_t, names = ['Direction','Frequency','Fluence','Time'])    
    
    data = np.zeros([len(direcs)*len(freqs)*len(fluences),len(time_arr)])

    backvs = []
    lowpow = min(fluences)
    
    re_backV = "^Background Voltage,-(\d+\.\d+)(.)V"
    re_backV2 = "^Background Voltage,-(\d+)(.)V"

    lp = 0
    for i, tup in enumerate(mi):
        direc, freq, fluence = tup
        fp = s_fps[direc,freq,fluence]
        if(fluence == lowpow):
            lp = load_trace(fp,50e-9).values
            data[i,:] = lp - lp
        else:
            d = load_trace(fp,50e-9).values

            data[i,:] = d - lp

        with open(fp) as p:
            for i, line in enumerate(p):
                if i == 11:
                    m = re.search(re_backV,line)
                    if m == None:
                        m = re.search(re_backV2,line)

                    if m.groups()[1] == 'm':
                        fac = 1e-3
                    elif m.groups()[1] == 'µ':
                        fac = 1e-6


                    backvs.append(float(m.groups()[0])*fac)

    data = data.flatten()
    
    s = pd.Series(data,index = mi_t)
    
    backvs = pd.Series(backvs,index = mi).xs(mi.levels[2][-1],level =2)
    
    return s,backvs


###Obsolete###

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
