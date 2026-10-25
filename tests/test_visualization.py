import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from fhn_control.models import FHNParameters, rhs
from fhn_control.visualization import plot_phase_portrait, animate_trajectory


def test_plot_phase_and_animation_html():
    p=FHNParameters(a=0.7,b=0.8,tau=12.5,I=0.9)
    t=np.linspace(0,40,50)
    sol=solve_ivp(lambda t,x:rhs(t,x,p),(0,40),[-1.5,0.1],t_eval=t)
    fig,ax=plt.subplots()
    plot_phase_portrait(ax,p,sol)
    fig.canvas.draw()
    plt.close(fig)
    anim=animate_trajectory(sol,p,max_frames=6)
    html=anim.to_jshtml(fps=10,default_mode='once')
    assert 'animation' in html.lower()
    plt.close(anim._fig)
