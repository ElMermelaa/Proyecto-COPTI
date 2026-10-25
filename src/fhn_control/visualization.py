"""Visualizaciones reutilizables para la dinámica de una neurona FHN.

Todas las funciones dibujan en memoria. No guardan archivos en disco.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from .models import equilibrium, rhs


def plot_phase_portrait(ax, p, sol=None, show_nullclines=True,
                        xlim=(-2.7, 2.7), ylim=(-1.8, 2.9),
                        tail_fraction=None):
    """Dibuja campo, nulclinas (opcional), equilibrio y trayectoria."""
    v=np.linspace(*xlim,21)
    w=np.linspace(*ylim,19)
    V,W=np.meshgrid(v,w)
    velocities=rhs(0,(V,W),p)
    speed=np.maximum(np.hypot(*velocities),1e-12)
    ax.quiver(V,W,velocities[0]/speed,velocities[1]/speed,alpha=.24)
    if show_nullclines:
        vv=np.linspace(*xlim,450)
        ax.plot(vv,vv-vv**3/3+p.I,lw=1.7,label=r'$\dot v=0$')
        if p.b>0:
            ax.plot(vv,(vv+p.a)/p.b,lw=1.7,label=r'$\dot w=0$')
        else:
            ax.axvline(-p.a,lw=1.7,label=r'$\dot w=0$')
    xstar=equilibrium(p)
    ax.scatter(*xstar,c='black',s=38,zorder=6,label='Equilibrio')
    if sol is not None:
        ax.plot(sol.y[0],sol.y[1],color='C2',lw=1.5,alpha=.9,label='Trayectoria')
        ax.scatter(*sol.y[:,0],color='C2',s=24,zorder=7,label='Inicial')
        if tail_fraction is not None and 0<tail_fraction<1:
            k=int((1-tail_fraction)*len(sol.t))
            ax.plot(sol.y[0,k:],sol.y[1,k:],color='C3',lw=2.6,label='Tramo final')
    ax.set(xlim=xlim,ylim=ylim,xlabel='$v$',ylabel='$w$')
    ax.grid(alpha=.15)
    ax.legend(fontsize=8,loc='best')
    return ax


def animate_trajectory(sol,p,max_frames=110,title='Evolución temporal'):
    """Devuelve una FuncAnimation con dos ejes y reproductor HTML interactivo.

    El llamador puede usar HTML(anim.to_jshtml(...)) dentro de Jupyter.
    No exporta video, GIF ni imágenes.
    """
    if not sol.success or len(sol.t)<2:
        raise ValueError('La integración debe ser exitosa y tener al menos dos tiempos')
    indices=np.unique(np.linspace(0,len(sol.t)-1,max_frames,dtype=int))
    fig,axs=plt.subplots(1,2,figsize=(9.8,4.1))
    phase,time_ax=axs
    phase.plot(sol.y[0],sol.y[1],color='C0',lw=1,alpha=.18)
    eq=equilibrium(p)
    phase.scatter(*eq,s=35,color='black',label='Equilibrio',zorder=5)
    moving_path,=phase.plot([],[],color='C0',lw=2.2,label='Recorrido')
    moving_dot,=phase.plot([],[],marker='o',ms=8,color='C3')
    vr=max(np.ptp(sol.y[0]),.5)
    wr=max(np.ptp(sol.y[1]),.5)
    phase.set(xlim=(min(sol.y[0].min(),eq[0])-.12*vr,
                     max(sol.y[0].max(),eq[0])+.12*vr),
              ylim=(min(sol.y[1].min(),eq[1])-.12*wr,
                    max(sol.y[1].max(),eq[1])+.12*wr),
              xlabel='$v$',ylabel='$w$',title='Plano de fases')
    phase.legend(loc='best',fontsize=8)
    time_ax.plot(sol.t,sol.y[0],alpha=.18,lw=1,color='C0')
    moving_v,=time_ax.plot([],[],color='C0',lw=2)
    marker,=time_ax.plot([],[],marker='o',ms=7,color='C3')
    time_ax.set(xlim=(sol.t[0],sol.t[-1]),
                ylim=(sol.y[0].min()-.12*vr,sol.y[0].max()+.12*vr),
                xlabel='t',ylabel='$v(t)$',title='Evolución de v(t)')
    for ax in axs: ax.grid(alpha=.2)
    fig.suptitle(title)
    fig.tight_layout()
    def update(frame):
        k=indices[frame]
        moving_path.set_data(sol.y[0,:k+1],sol.y[1,:k+1])
        moving_dot.set_data([sol.y[0,k]],[sol.y[1,k]])
        moving_v.set_data(sol.t[:k+1],sol.y[0,:k+1])
        marker.set_data([sol.t[k]],[sol.y[0,k]])
        return moving_path,moving_dot,moving_v,marker
    return FuncAnimation(fig,update,frames=len(indices),interval=65,blit=False,repeat=False)
