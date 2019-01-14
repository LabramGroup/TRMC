import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, FuncFormatter


class exp_formatter(): 
    """used to format exponentials of ticks"""
    def __init__(self,exponent):
        self.exponent = exponent
    
    def format_func(self,value, tick_number):
        return ("{:0=1.0f}").format(value/10**self.exponent)
#         return ("{:0=1.0f}e" + str(self.exponent)).format(value/10**self.exponent)



def dvcolorplot(sweep, dvs , levels = list(np.arange(-3.1e-3,3.1e-3,1e-5))):
    """
    Color plot of delta v vs time and freq
    
    Data input  is one cavity sweep and a  multindex Series for deltaVs with time and freq as levels
    """
    
    freq = dvs.index.levels[1]
    time = dvs.index.levels[0]
    z = dvs.values.reshape(len(time),len(freq))

    xi, yi = np.meshgrid(freq,time)

    fig , axes = plt.subplots(2,1 , sharex = True,constrained_layout=True)

    cs = axes[0].contourf(xi,yi,z, levels, cmap='seismic')
    cb = fig.colorbar(cs, ax = axes[0])


    expf = exp_formatter(-9)
    axes[0].yaxis.set_major_formatter(FuncFormatter(expf.format_func))
    axes[0].set_ylabel('Time (ns)')
    axes[1].set_xlabel('Freq (Hz)')
    axes[1].set_ylabel('Normalized\n Reflectivity')
    cb.set_label('Voltage (V)')

    axes[1].plot(sweep)

    return fig, axes    


def absplot(dvs):
    """
    Plots traces with negative integral as positive but red color

    Data is multindex with time and freq as levels
    """

    fig, ax = plt.subplots()

    for freq in dvs.index.levels[1]:
        trace = dvs[:,freq]
        if len(trace.shape) == 1:
            if np.trapz(trace) > 0:
                color = 'b'
                zorder = 2
            else:
                color = 'r'
                zorder = 1
            ax.plot(abs(trace), color = color, zorder = zorder )


    ax.set_yscale('log')
    ax.set_ylim([1e-7,5e-3])
    ax.set_ylabel('Voltage (V)')
    ax.set_xlabel('Time (ns)')

    expf = exp_formatter(-9)
    ax.xaxis.set_major_formatter(FuncFormatter(expf.format_func))
    return fig, ax

def vsplotxr(timesel, dvs, vss = None, fits = None, v0 = None, v0_fit = None):
    timesel = timesel *1e-9
    if fits is not None:
        times = fits.indexes['time']
        sample = fits.sample.values
    elif vss is not None:
        times = vss.indexes['time']
        sample = vss.sample.values
    
    timesel = min(times, key=lambda x:abs(x-timesel))

    fig, axes = plt.subplots(2,1, figsize = (7,10), sharex = True)
    axes[0].axhline(0, color = 'gray', linestyle = '--')
    axes[0].plot(dvs.loc[:,timesel].to_series(), marker = 'o')

    if v0 is not None:
        v0 = v0.to_series()
        axes[1].plot(v0, label = 'cavity sweep', marker = 'o')
        axes[1].axvline(v0.idxmin(), color = 'gray', linestyle = '--')
        axes[0].axvline(v0.idxmin(), color = 'gray', linestyle = '--')
    if v0_fit is not None:
        axes[1].plot(v0_fit.to_series(), label = 'cavity sweep fit')
        
    
    if vss is not None:
        axes[1].plot(vss.loc[:,timesel].to_series(), color = 'r', label = 'cavity sweep + deltaV' ,  marker = 'o')
#         print(vss.loc[:,timesel].to_series())
    if fits is not None:
        axes[1].plot(fits.loc[:,timesel].to_series(), label = 'cavity sweep + deltaV fit')
#         print(fits.loc[:,timesel].to_series())

    axes[0].set_ylabel('Delta V (V)')
    axes[1].set_ylabel('$V_r$ (V)')
    axes[1].set_xlabel('Frequency (Hz)')

    # axes[1].set_xlim([8.525e9,8.545e9])
    # axes[1].set_ylim([0.4,0.8])

    fig.suptitle('Delta V taken at ' + str(timesel) + 's for sample ' + str(sample))
    
    return fig, axes
