import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, FuncFormatter

from matplotlib import animation, rc
from IPython.display import HTML


class exp_formatter(): 
    """used to format exponentials of ticks"""
    def __init__(self,exponent):
        self.exponent = exponent
    
    def func(self,value, tick_number):
        return ("{:0=1.0f}").format(value/10**self.exponent)
#         return ("{:0=1.0f}e" + str(self.exponent)).format(value/10**self.exponent)



def dvcolorplot(sweep, dvs , levels = list(np.arange(-3.1e-3,3.1e-3,1e-5))):
    """
    Color plot of delta v vs time and freq
    
    Data input  is one $V_{bg}(\omega)$ and a  multindex Series for deltaVs with time and freq as levels
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

    fig, axes = plt.subplots(2,1, figsize = (10,15), sharex = True)

    axes[0].axhline(0, color = 'gray', linestyle = '--')
    axes[0].plot(dvs.loc[:,timesel].to_series(), marker = 'o', label = '$\Delta V(\omega)$')

    if v0 is not None:
        v0 = v0.to_series()
        axes[1].plot(v0, label = '$V_{bg}(\omega)$', marker = 'o')
        axes[1].axvline(v0.idxmin(), color = 'gray', linestyle = '--')
        axes[0].axvline(v0.idxmin(), color = 'gray', linestyle = '--')
    if v0_fit is not None:
        axes[1].plot(v0_fit.to_series(), label = '$V_{bg}(\omega)$ fit')
        
    
    if vss is not None:
        axes[1].plot(vss.loc[:,timesel].to_series(), color = 'r', label = '$V_{bg}(\omega)$ +  $\Delta V(\omega)$' ,  marker = 'o')
#         print(vss.loc[:,timesel].to_series())
    if fits is not None:
        axes[1].plot(fits.loc[:,timesel].to_series(), label = '$V_{bg}(\omega)$ + $\Delta V(\omega)$ fit')
#         print(fits.loc[:,timesel].to_series())

    axes[0].set_ylabel('$\Delta V(\omega)$ (V)')
    axes[1].set_ylabel('$V_r$ (V)')
    axes[1].set_xlabel('Frequency (Hz)')

    # axes[1].set_xlim([8.525e9,8.545e9])
    # axes[1].set_ylim([0.4,0.8])

    fig.suptitle('$\Delta V(\omega)$ taken at ' + str(timesel) + 's for sample ' + str(sample))
    
    return fig, axes


# First set up the figure, the axis, and the plot element we want to animate

def sweepfitanim(dst,interval = 50):

    fittimes = dst.indexes['time']
    RawData_Frames = dst['vss'].loc[:,fittimes]
    RawData_Frames_fit = dst['fits'].loc[:,fittimes]

    # interval = 50
    xs = dst.indexes['freq']

    fig = plt.figure()
    ax = plt.axes(xlim=(xs[0], xs[-1]), ylim = (0.005,0.025))

    line, = ax.plot([], [], lw=2, marker = 'o')
    line_fit, = ax.plot([], [], lw=2, color = 'r')

    time_template = 'Time = %.1fns'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    fig.tight_layout()

    # initialization function: plot the background of each frame
    def init():
        line.set_data([], [])
        line_fit.set_data([], [])
        time_text.set_text('')
        ax.set_ylabel("Voltage (V)")
        ax.set_xlabel("Frequency (Hz)") 
        return line,line_fit

    # animation function.  This is called sequentially
    def animate(i):
        line.set_data(xs, RawData_Frames.loc[:,fittimes[i]])
        line_fit.set_data(xs, RawData_Frames_fit.loc[:,fittimes[i]])
        time_text.set_text(time_template % int(fittimes[i]*1e9))
        return line, time_text

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(fittimes), interval = interval, blit=True)

    rc('animation', html='html5')

    return anim

def redbluetransient(ax,data,f0):
    freqs = data.indexes['freq']
    labeledblue = False
    labeledred = False
    for freq in freqs:
        if freq < f0:
            color = 'b'
            if labeledblue:
                label = None
            else:
                label = 'Frequency Below Resonance'
                labeledblue = True
        elif freq > f0:
            color = 'r'
            if labeledred:
                label = None
            else:
                label = 'Frequency Above Resonance'
                labeledred = True
        elif freq == f0:
            color = 'g'
            label = 'On resonance'
        trace = data.sel(freq = freq)
        ax.plot(trace.to_series(), color = color, label = label)

    ax.xaxis.set_major_formatter(FuncFormatter(exp_formatter(-9).func))

    # ax.set_xlim(200e-9,400e-9)
    # ax.set_ylim([-0.002,0.002])
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Voltage (V)')


    return ax