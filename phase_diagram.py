import numpy as np
from phase_portrait import roots_fun

def det_phase (G,C,D,B):
    roots = roots_fun(C,G,D,B)
    if len(roots)==3 and max(roots)>=C/2:
        return 'bistability'
    else:
        root = min(roots)
        A = root
        R = D/(D-B*A)
        Jac = [[-D+B*A,B*R],
               [-A,-R+2*A*(1-A/C)*G-G*A**2/C]]
        w,v = np.linalg.eig(Jac)
        if root<1/np.sqrt(G) and R>G*C/4:
            return 'single healthy'
        elif root<1/np.sqrt(G):
            return 'flare ups'
        elif (np.iscomplex(w[0]) and np.real(w[0])>0) or (w[0]>0 and w[1]>0):
            return 'oscillations'
        elif (np.iscomplex(w[0]) and np.real(w[0]<0)) or np.iscomplex(w[0])==False:
            return 'single sick'
        else:
            return ''

def compute_phase_diagram (G_min,G_max,B_min,B_max,C,D):

    Gs = np.logspace(np.log10(G_min),np.log10(G_max),500)
    Bs = np.logspace(np.log10(B_min), np.log10(B_max),500)
    phases = {}

    for G in Gs:
        for B in Bs:
            phase = det_phase(G,C,D,B)
            if phase!='':
                phases[phase].append((B,G))
    
    return phases

def create_phase_diagram_figure (phases,colors,G_min,G_max,B_min,B_max,phases_names):
    import matplotlib.pyplot as plt
    fig,axs = plt.subplots(figsize = (10,10))
    for phase,color in zip(phases_names,colors):
        axs.scatter(*zip(*phases[phase]),label = phase,color=color,s = 1, alpha = 1)
        
    axs.set_xscale('log')
    axs.set_yscale('log')
    axs.set_xlim(B_min,B_max)
    axs.set_ylim(G_min,G_max)
    axs.set_xlabel('B',fontsize = 20)
    axs.set_ylabel('G',fontsize = 20)
    axs.grid(False)
    #axs.axis('square')
    axs.legend(loc = 'upper right')
    axs.tick_params(axis='both', which='major', labelsize=12)
    axs.tick_params(axis='both', which='minor', labelsize=10)
    axs.spines[['right', 'top']].set_visible(False)
    return fig,axs

if __name__=="__main__":
    phases_names = ['bistability','single healthy','flare ups','oscillations','single sick']
    colors_dict =  {'bistability':'#fdae61','single healthy':'#2c7bb6','flare ups':'#abd9e9','oscillations':'#ffffbf','single sick':'#d7191c'}           
    colors = ['#fdae61','#2c7bb6','#abd9e9','#ffffbf','#d7191c']
    G_min = 1e-3
    G_max = 1
    B_min = 1e-4
    B_max = 1
    C=1000
    D = 0.15
    phases = compute_phase_diagram(G_min,G_max,B_min,B_max,C,D)
    create_phase_diagram_figure(phases,colors,G_min,G_max,B_min,B_max,phases_names)