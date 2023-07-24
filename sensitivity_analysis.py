import numpy as np
import matplotlib.pyplot as plt
from mpl_default_settings import load_def_settings

def compute_sensitivity_analysis (epsilons,q,r_baseline,std,funs):
    q_m = q-std
    q_p = q+std
    fold_change = [1+e for e in epsilons]
    ys = []
    y_ms = []
    y_ps = []
    for i,fun in enumerate(funs):
        ys.append([fun(_,r_baseline,q) for _ in epsilons])
        y_ms.append([fun(_,r_baseline,q_m) for _ in epsilons])
        y_ps.append([fun(_,r_baseline,q_p) for _ in epsilons])
    
    return fold_change, ys,y_ms,y_ps

def create_figure_sensitivity_analysis (fold_change,ys,y_ms,y_ps,lbls,colors):
    load_def_settings()
    fig,ax = plt.subplots()

    for i,y in enumerate(ys):
        ax.plot(fold_change,y, label = lbls[i],color = colors[i],zorder = len(funs)-i)
        ax.fill_between(fold_change,y_ms[i],y_ps[i],color = colors[i],alpha = 0.2,zorder = len(funs)-i)
            
    ax.legend(bbox_to_anchor=(1.1, 1.05),frameon=True)
    ax.set_xlabel('fold change in parameter')
    ax.set_ylabel('$r/r_{baseline}$')
    ax.set_yscale('log')
    ax.spines[['top','right']].set_visible(False)
    ax.set_xlim(min(fold_change),max(fold_change))
    return fig,ax

if __name__=="__main__":
    change_h_lR = lambda eps,r_baseline,q:r_baseline**(3*eps)*np.exp(-2*eps*q)
    change_mR = lambda eps,r_baseline,q:r_baseline**(-3*eps)*np.exp(2*eps*q)
    change_lA = lambda eps,r_baseline,q: np.exp(q*eps)
    change_gamma = lambda eps,r_baseline,q:r_baseline**(-2*eps)*np.exp(3*eps*q)
    change_sigma = lambda eps,r_baseline,q:r_baseline**(eps)*np.exp(-2*eps*q)
    change_beta_C_mA = lambda eps,r_baseline,q:1
    funs = [change_gamma,change_h_lR,change_mR,change_sigma,change_lA,change_beta_C_mA]

    epsilons = np.linspace(-0.15,0.15,100)
    q = 10
    r_baseline = 0.64
    std = 4
    lbls = ['$\gamma$','$h,l_R$','$m_R$','$\sigma$','$l_A$','$\\beta,C,m_A$']
    colors = ['#313695','#4575b4','#4575b4','#74add1','#abd9e9','#e0f3f8']
    fold_change,ys,y_ms,y_ps = compute_sensitivity_analysis(epsilons, q, r_baseline, std, funs)
    create_figure_sensitivity_analysis(fold_change, ys, y_ms, y_ps, lbls, colors)