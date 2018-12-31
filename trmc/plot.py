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


