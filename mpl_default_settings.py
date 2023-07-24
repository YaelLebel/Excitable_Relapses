import matplotlib.pyplot as plt
import matplotlib as mpl

def load_def_settings ():
    plt.rcParams["figure.figsize"] = (2, 2)
    plt.rcParams['lines.linewidth'] = 2
    plt.rc('xtick', labelsize=7) 
    plt.rc('ytick', labelsize=7) 
    plt.rc('font', size=7) 