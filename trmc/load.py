import os
import pandas as pd
import numpy as np
import re
import xarray as xr
from IPython.display import clear_output


def loadsweep(fp,defaultV = 0.025):
    """Load in cavity sweep. defaultV if bad csv file saved. """ 
    df = pd.read_csv(fp,index_col=False)
    
    if 'Experimental R' in df.columns:
        df = pd.read_csv(fp,index_col=False)
        df = df.set_index(df['f(Ghz)'])
        s = df[' Vsignal(V)']
    else:
        df = pd.read_csv(fp, index_col = 0, skiprows = 4)
        s = df['Reflectivity']*defaultV
    return s



#Load in cavity sweeeps
def sweeps2ds(fps):
    """load in all cavity sweeps in filepath dict"""
    regex = 'Sweep_(\d+)ms(.+)exp.csv'

    das = []
    for samp in fps:
        direc = fps[samp]
        fns = os.listdir(direc)
        for fn in fns:
            m = re.search(regex,fn)
            if m is None:
                pass
            else:
                fp = os.path.join(direc,fn)
                swtime = int(m.groups()[0])
                tc = m.groups()[1]            

                s = loadsweep(fp)
                s = s.rename(s.name.replace(' ', ''))
                s.index = s.index.rename('freq')
                da = xr.DataArray.from_series(s)
                da = da.assign_coords(sample = samp).expand_dims('sample')
                da = da.assign_coords(tc = tc).expand_dims('tc')
                da = da.assign_coords(swtime= swtime).expand_dims('swtime')
                das.append(da)

    ds = xr.merge(das)
    return ds

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
    file_re = '.*Filter=\d+_Fluence=(.+?)_data.csv'

    flist = []

    for folder in folders:
        m_folder = re.search(folder_re,folder)
        freq = float(m_folder.groups(0)[0])*1e9
        direction = m_folder.groups(0)[1]
        folderpath = os.path.join(direc,folder)
        fns = os.listdir(folderpath)
        for fn in fns:
            #if file[0] == 'F':
            m_file = re.search(file_re, fn)
            if m_file is None:
                clear_output()
                print("no match for file " + fn)
            else:
                fluence = float(m_file.groups(0)[0])

                fp = os.path.join(folderpath,fn)

                miarray.append((direction,freq,fluence))
                flist.append(fp)


    mi = pd.MultiIndex.from_tuples(miarray, names = ['direction','freq','fluence'])

    s_fps = pd.Series(flist, index = mi)
    
    return s_fps

def freqdcs_flist(direc):
    """Creates a multindexed Series of filepaths from a frequency dark cavity sweep folder"""
    folders = os.listdir(direc)
    miarray = []
    folder_re = '^(\d+\.\d+)GHz_(.+?)'
    # file_re = '.*Filter=\d+_Fluence=(.+?)_data.csv'

    flist = []

    for folder in folders:
        m_folder = re.search(folder_re,folder)
        freq = float(m_folder.groups(0)[0])*1e9
        direction = m_folder.groups(0)[1]
        folderpath = os.path.join(direc,folder)
        fns = os.listdir(folderpath)
        for fn in fns:
            if fn == 'FreqSweep_DarkCavitySweep_exp.csv':
                fp = os.path.join(folderpath,fn)
                miarray.append((direction,freq))
                flist.append(fp)


    mi = pd.MultiIndex.from_tuples(miarray, names = ['direction','freq'])

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

    mi_t = pd.MultiIndex.from_tuples(miarray_t, names = ['direction','freq','fluence','time'])    
    
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
            try:
                data[i,:] = d - lp
            except:
                print('subtraction failed for ' + fp)

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
