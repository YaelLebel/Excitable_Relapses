import numpy as np
from model_equations import model,Adot,Rdot,A_nullcline,R_nullcline
import matplotlib.pyplot as plt

def roots_fun (C,G,D,B):
    roots = np.roots([B*G/C,-B*G-D*G/C,D*G,-B-D,D])
    roots = roots[np.iscomplex(roots)==False]
    roots = roots[roots>0]
    roots = roots[roots<=D/B]
    roots = np.sort(roots)
    return roots

def create_phase_portrait(C,G,D,B):  
    #setting up nullclines and vector field
    xmax = C
    x = np.linspace(-C,xmax,1000)
    ymax = 1.1*G*C/4
    y = np.linspace(-ymax,ymax,1000)
    X,Y = np.meshgrid(x,y)
    u = Adot(X,Y,G,C,D,B)
    v = Rdot(X,Y,G,C,D,B)
    X_A = np.linspace(1e-5,C,10000)
    A_nc = A_nullcline(X_A,C,G,D,B)
    X_R = np.linspace(-10,D/B,10000)
    R_nc = R_nullcline(X_R,C,G,D,B)

    from mpl_default_settings import load_def_settings
    
    load_def_settings()

    fig,ax = plt.subplots()
    
    ax.streamplot(X,Y,u,v, density = 0.5,broken_streamlines=False, color = '#C5C9C7',zorder = -1)
    
    ax.plot(X_A,A_nc,label = 'A nullcline',color = '#2b83ba',zorder = -1)
    ax.plot(X_R,R_nc, label = 'R nullcline',color = '#A60628',zorder = -1)
    ax.set_xlim(-10,xmax)
    ax.set_ylim(0,ymax)
    ax.set_xlabel('$A / A0$')
    ax.set_ylabel('$R / R0$')
    roots = roots_fun(C,G,D,B)
    
    ax.scatter(roots[0],D/(D-B*roots[0]),s = 50,zorder=1,color = 'black')
    ax.scatter(roots[1],D/(D-B*roots[1]),facecolors='none' ,s = 50, color = 'black',zorder = 1)
    ax.scatter(roots[2],D/(D-B*roots[2]),facecolors='none',s = 50, color = 'black',zorder = 1)
    return fig,ax

def create_inset_fixed_points (C,G,D,B):
    fig,ax = plt.subplots(figsize = (1,1))
    xmin = 1e-1
    xmax = 29
    x1 = np.linspace(xmin,xmax,1000)
    A_nullcline_inset = [(1+G*(x**2)*(1-(x/C)))/(x) for x in x1]
    R_nullcline_inset = [D/(D-B*x) for x in x1]
    ax.plot(x1,R_nullcline_inset, label = 'R nullcline',color = '#A60628',zorder = -1)
    ax.plot(x1,A_nullcline_inset,label = 'A nullcline',color = '#2b83ba',zorder = -1)
    roots = roots_fun(C,G,D,B)
    ax.scatter(roots[0],D/(D-B*roots[0]),s = 50,zorder=1,color = 'black')
    ax.scatter(roots[1],D/(D-B*roots[1]),facecolors='white' ,s = 50, color = 'black',zorder = 1)
    ax.scatter(roots[2],D/(D-B*roots[2]),facecolors='white',s = 50, color = 'black',zorder = 1)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(5e-1,1e2)
    ax.set_xlim(xmin,xmax+10)
    ax.spines[['right','top']].set_visible(False)


def create_phase_portrait_with_trajectory (C,G,D,B,Tfinal,num_steps):
    from scipy.integrate import odeint

    t = np.linspace(0,Tfinal,num_steps)
    R_icc = 1
    Ath = 1/G
    A_icc = 2*Ath
    sol = odeint(model,[A_icc,R_icc],t,args = (G,C,D,B))
    fig,ax = create_phase_portrait(C, G, D, B)
    line = ax.plot(sol[:,0],sol[:,1], linestyle = '--',color = 'black')[0]
    pos = np.arange(0,int(num_steps*0.2),75)
    for pos in pos:
        line.axes.annotate('',
            xytext=(sol[:,0][pos], sol[:,1][pos]),
            xy=(sol[:,0][pos+1], sol[:,1][pos+1]),
            arrowprops=dict(arrowstyle="-|>", color=line.get_color()),
            size=15)
    legend_properties = {'weight':'bold'}
    ax.legend(bbox_to_anchor=(1.1, 1.05),prop = legend_properties,frameon=True)

if __name__=="__main__":
    Tfinal = 100
    num_steps = 10000
    
    C = 1000
    G = 0.25
    D = 0.15
    B = 0.001
    
    create_phase_portrait(C, G, D, B)
    create_inset_fixed_points(C, G, D, B)
    create_phase_portrait_with_trajectory(C, G, D, B, Tfinal, num_steps)