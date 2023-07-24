import numpy as np
import matplotlib.pyplot as plt
from mpl_default_settings import load_def_settings
from rfp_hazard_bootstrap import rfp_hazard

def rfp_hazard_simulation (filename):
    diffs_all = np.loadtxt(filename)
    rfp_all,hazard_all,bins = rfp_hazard(diffs_all,bins=50)

    return rfp_all, hazard_all, bins

def create_figure_rfp_hazard_simulations (rfp_all,hazard_all,bins):
    load_def_settings()
    color_s = '#fc8d59'
    color_h = '#91bfdb'
    fig,ax = plt.subplots()
    ax.plot(bins[1:], rfp_all, color = color_s,label="survival")
    ax.plot(bins[2:-1], hazard_all, color = color_h,label="hazard")
    ax.legend()
    ax.spines[['top','right']].set_visible(False)
    ax.set_xlabel('Time [days]')
    return fig,ax

if __name__ =='__main__':
    filename = "diffsG01D015B0001.txt"
    rfp_all,hazard_all,bins = rfp_hazard_simulation (filename)
    create_figure_rfp_hazard_simulations(rfp_all,hazard_all,bins)
