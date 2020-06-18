#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 15:16:07 2020

@author: Francesco Di Lauro
@mail: F.Di-Lauro@sussex.ac.uk
Copyright 2020 Francesco Di Lauro. All Rights Reserved.
See LICENSE file for details
"""

from scipy.integrate import  odeint
import numpy as np
##############

import matplotlib.pyplot as plt
from matplotlib import rc
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def sir(v,t, beta, gamma, int_time, duration, c):
    '''
         ODE for SIR with control:
             The system to solve is:
                 ds/dt =  -beta i s
                 di/dt = beta i s - gamma i
                 dr/dt = gamma i
             Parameters: v: VECTOR, current value of the dependent variables (v[0]=S, v[1]=I, v[2]=R).
                         t: float, indipendent variable (time).
                         beta: float, infection parameter.
                         gamma: float, recovery parameter.
                         int_time: float,
                         c_red: float, 0<c_red<1 PERCENTAGE of reduction, c_red = 0 means no reduction occurs, otherwise reduction of tau by (1-c)*tau
                         c_duration: float, duration of control measures.
             Returns: value of vdot
    '''
    vdot = np.empty_like(v)
    #If you use vdot = v, vdot is not a copy of v. The two names now refer to the same object.
    #So when you start modifying vdot in the function Diffeq, you are actually modifying the input argument.
    #Apparently that affects the behavior of odeint.
    if t> int_time and t < int_time+duration:
        beta_i = (1-c)*beta
    else:
        beta_i = beta

    vdot[0] = -beta_i*v[0]*v[1]  #susceptible
    vdot[1] = beta_i*v[0]*v[1] - gamma*v[1]  #infected
    vdot[2] = gamma*v[1]  #recovered

    return vdot

def sir_odeint( beta, gamma, int_time, duration, c, tauf,nt, incond):

        t = np.linspace(0,tauf,nt)
        #args=(ak,dk)
        y= odeint(sir, incond,t, args=(beta, gamma, int_time, duration, c))
        return y
def SversusRfig():  #I renamed this from fig1 --- Joel
    fig,ax = plt.subplots(1,3, figsize=(10,2.8),sharey=True)
    fig.subplots_adjust(wspace=0.3)

    ax = ax.ravel()
    betav=[0.5,2.0,4.0]
    gamma=1
    Sv = np.linspace(.08,0.9,5)
    
    tauf= 20
    nt = 8000
    t = np.linspace(0,tauf,nt)
    hwidth = 0.05
    size=10
    label=[r'$(a)$',r'$(b)$',r'$(c)$']
    for i, sub in enumerate(ax):
        beta = betav[i]
        for S in Sv:
            incond = np.array([S,1-S,0])
    
            y_noint = sir_odeint(beta,gamma,-2,-2,0, tauf, nt, incond)
            sub.plot(y_noint[:,2],y_noint[:,0],color='gray', linestyle='--')      
            if i>0:
                sub.scatter(1-1.0/beta, 1.0/beta, color='r')
            #plt.scatter(0,S)
            adx0= np.argwhere( y_noint[:,2] > 0.15)[0][0]
            
            sub.text(0.9, 0.9,label[i], size = 12)

            arrow0 = y_noint[:,2][adx0+1], y_noint[:,0][adx0+1], -y_noint[:,2][adx0]+y_noint[:,2][adx0+1], -y_noint[:,0][adx0]+y_noint[:,0][adx0+1]
    
    
            sub.arrow(*arrow0, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k') 
            try:
                adx1= np.argwhere(y_noint[:,2] > 0.3)[0][0]
                arrow1 = y_noint[:,2][adx1+1], y_noint[:,0][adx1+1], -y_noint[:,2][adx1]+y_noint[:,2][adx1+1], -y_noint[:,0][adx1]+y_noint[:,0][adx1+1]
                sub.arrow(*arrow1, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k')  
            except:
                pass
            try:
                adx2 =  np.argwhere( y_noint[:,2] > 0.45)[0][0]
                arrow2 = y_noint[:,2][adx2+1], y_noint[:,0][adx2+1], -y_noint[:,2][adx2]+y_noint[:,2][adx2+1], -y_noint[:,0][adx2]+y_noint[:,0][adx2+1]
                sub.arrow(*arrow2, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k')
            except:
                pass
                
            try:
                adx3 = np.argwhere( y_noint[:,2] > 0.6)[0][0]
                arrow3 = y_noint[:,2][adx3+1], y_noint[:,0][adx3+1], -y_noint[:,2][adx3]+y_noint[:,2][adx3+1], -y_noint[:,0][adx3]+y_noint[:,0][adx3+1]
                sub.arrow(*arrow3, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k')
            except:
                pass
            try:
                adx5 = np.argwhere( y_noint[:,2] > 0.9)[0][0]
                arrow5 = y_noint[:,2][adx5+1], y_noint[:,0][adx5+1], -y_noint[:,2][adx5]+y_noint[:,2][adx5+1], -y_noint[:,0][adx5]+y_noint[:,0][adx5+1]
                sub.arrow(*arrow5, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k')            
            except:
                pass
            try:    
                adx4 = np.argwhere( y_noint[:,2] > 0.73)[0][0]
                arrow4 = y_noint[:,2][adx4+1], y_noint[:,0][adx4+1], -y_noint[:,2][adx4]+y_noint[:,2][adx4+1], -y_noint[:,0][adx4]+y_noint[:,0][adx4+1]
                sub.arrow(*arrow4, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k') 
            except:
                pass
        
        if i==1:
            incond = np.array([0.81,0.01,0.18])
    
            y_noint = sir_odeint(beta,gamma,-2,-2,0, tauf, nt, incond)
            sub.plot(y_noint[:,2],y_noint[:,0],color='gray', linestyle='--')    

            try:
                adx1= np.argwhere(y_noint[:,2] > 0.3)[0][0]
                arrow1 = y_noint[:,2][adx1+1], y_noint[:,0][adx1+1], -y_noint[:,2][adx1]+y_noint[:,2][adx1+1], -y_noint[:,0][adx1]+y_noint[:,0][adx1+1]
                sub.arrow(*arrow1, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k')  
            except:
                pass
            try:
                adx2 =  np.argwhere( y_noint[:,2] > 0.45)[0][0]
                arrow2 = y_noint[:,2][adx2+1], y_noint[:,0][adx2+1], -y_noint[:,2][adx2]+y_noint[:,2][adx2+1], -y_noint[:,0][adx2]+y_noint[:,0][adx2+1]
                sub.arrow(*arrow2, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k')
            except:
                pass
                
            try:
                adx3 = np.argwhere( y_noint[:,2] > 0.6)[0][0]
                arrow3 = y_noint[:,2][adx3+1], y_noint[:,0][adx3+1], -y_noint[:,2][adx3]+y_noint[:,2][adx3+1], -y_noint[:,0][adx3]+y_noint[:,0][adx3+1]
                sub.arrow(*arrow3, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k')
            except:
                pass
            try:
                adx5 = np.argwhere( y_noint[:,2] > 0.9)[0][0]
                arrow5 = y_noint[:,2][adx5+1], y_noint[:,0][adx5+1], -y_noint[:,2][adx5]+y_noint[:,2][adx5+1], -y_noint[:,0][adx5]+y_noint[:,0][adx5+1]
                sub.arrow(*arrow5, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k')            
            except:
                pass
            try:    
                adx4 = np.argwhere( y_noint[:,2] > 0.73)[0][0]
                arrow4 = y_noint[:,2][adx4+1], y_noint[:,0][adx4+1], -y_noint[:,2][adx4]+y_noint[:,2][adx4+1], -y_noint[:,0][adx4]+y_noint[:,0][adx4+1]
                sub.arrow(*arrow4, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k') 
            except:
                pass            
        x=np.linspace(0,1,100)
        sub.plot(x,1-x, linestyle='--', color='k')
    
        sub.set_xlim(0,1)
        sub.set_ylim(0,1)
        l2 = np.array([0.45,0.6])

        trans_angle = plt.gca().transData.transform_angles(np.array((-45,)),
                                                   l2.reshape((1, 2)))[0]
        sub.text(l2[0],l2[1], r'$S+R=1$', fontsize=size,rotation=trans_angle, rotation_mode='anchor')
        
        sub.set_title(r"$\mathcal{R}_0 = %.1f$"%beta, fontsize=size, position=(0.5, 0.86))
  
        sub.tick_params(labelsize=8, which='major')
        sub.tick_params(labelsize=2, which='minor')

        sub.set_xticks([0,0.3,0.6,0.9])
        sub.set_xticklabels([r'$0$',r'$0.3$',r'$0.6$',r'$0.9$'])
        sub.set_yticks([0,0.2,0.4,0.6,0.8,1])
        sub.set_yticklabels([r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$',r'$1$'])
        if i==0:
            sub.set_ylabel(r"$S$", size=11)
        if i==1:
            sub.set_xlabel(r"$R$", labelpad=-1,size=11)
     
    
    incondv = [[00.85,0.05,0.1],[0.77,0.03,0.2], [1-0.03-0.35, 0.03, 0.35]]
    for i,incond in enumerate(incondv):
        y_noint = sir_odeint(beta,gamma,-2,-2,0, tauf, nt, incond)
        ax[2].plot(y_noint[:,2],y_noint[:,0],color='gray', linestyle='--')      
        if i<2:
            try:
                adx1= np.argwhere(y_noint[:,2] > 0.3)[0][0]
                arrow1 = y_noint[:,2][adx1+1], y_noint[:,0][adx1+1], -y_noint[:,2][adx1]+y_noint[:,2][adx1+1], -y_noint[:,0][adx1]+y_noint[:,0][adx1+1]
                ax[2].arrow(*arrow1, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k')  
            except:
                pass
        try:
            adx2 =  np.argwhere( y_noint[:,2] > 0.45)[0][0]
            arrow2 = y_noint[:,2][adx2+1], y_noint[:,0][adx2+1], -y_noint[:,2][adx2]+y_noint[:,2][adx2+1], -y_noint[:,0][adx2]+y_noint[:,0][adx2+1]
            ax[2].arrow(*arrow2, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k')
        except:
            pass
            
        try:
            adx3 = np.argwhere( y_noint[:,2] > 0.6)[0][0]
            arrow3 = y_noint[:,2][adx3+1], y_noint[:,0][adx3+1], -y_noint[:,2][adx3]+y_noint[:,2][adx3+1], -y_noint[:,0][adx3]+y_noint[:,0][adx3+1]
            ax[2].arrow(*arrow3, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k')
        except:
            pass
        try:
            adx5 = np.argwhere( y_noint[:,2] > 0.9)[0][0]
            arrow5 = y_noint[:,2][adx5+1], y_noint[:,0][adx5+1], -y_noint[:,2][adx5]+y_noint[:,2][adx5+1], -y_noint[:,0][adx5]+y_noint[:,0][adx5+1]
            ax[2].arrow(*arrow5, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k')            
        except:
            pass
        try:    
            adx4 = np.argwhere( y_noint[:,2] > 0.73)[0][0]
            arrow4 = y_noint[:,2][adx4+1], y_noint[:,0][adx4+1], -y_noint[:,2][adx4]+y_noint[:,2][adx4+1], -y_noint[:,0][adx4]+y_noint[:,0][adx4+1]
            ax[2].arrow(*arrow4, shape='full', lw=0, length_includes_head=True, head_width=hwidth, color='k') 
        except:
            pass    
    plt.savefig("SversusRfig.eps",forma='eps', dpi=500)                
    

if __name__ == '__main__':
    
 
    SversusRfig()
    '''
    another interesting figure
    
    beta=2
    gamma=1.0
    c = 0.75
    tauf=16
    int_time = 3
    duration = 1
    nt = 5000
    t = np.linspace(0,tauf,nt)

    incond = np.array([0.99,0.01,0])
    
    y_noint = sir_odeint(beta,gamma,0,0,0, tauf, nt, incond)
    plt.figure()
    plt.plot(t,y_noint[:,1], label = 'no intervention', color='k')
    y_1 = sir_odeint(beta,gamma,1,duration,c, tauf, nt, incond)
    plt.plot(t,y_1[:,1], label = r'$t=1$', color='blue')
    y_3 =sir_odeint(beta,gamma,3,duration,c, tauf, nt, incond)
    plt.plot(t,y_3[:,1], label = r'$t=3$', color='brown')
    y_5 = sir_odeint(beta,gamma,4,duration,c, tauf, nt, incond)
    plt.plot(t,y_5[:,1], label = r'$t=4$', color='red')

    y_8 = sir_odeint(beta,gamma,8,duration,c, tauf, nt, incond)
    plt.plot(t,y_8[:,1], label = r'$t=8$', color='green')

    plt.xlabel(r"$t$")
    plt.ylabel(r"$I(t)$")
    plt.legend()    
    
    plt.figure()

    alpha=1
    plt.plot(t,y_1[:,0], label = r'$t=1$', color='blue')
    plt.plot(t,y_1[:,2], alpha=alpha,linestyle ='--',color='blue')


    plt.vlines(1, 0.85, 1, linestyle ='-', color='blue')
    plt.vlines(1+duration, 0.85, 1,linestyle ='--', color='blue')
    plt.plot(t,y_3[:,0], label = r'$t=3$', color='brown')
    plt.plot(t,y_3[:,2], alpha=alpha,linestyle ='--', color='brown')

    plt.vlines(3, 0.65, 0.8, linestyle ='-', color='brown')
    plt.vlines(3+duration, 0.65, 0.8,linestyle ='--', color='brown') 
    
    plt.plot(t,y_5[:,0], label = r'$t=4$', color='r')
    plt.plot(t,y_5[:,2],alpha=alpha,linestyle ='--', color='r')

    plt.vlines(4, 0.45, 0.66, linestyle ='-', color='r')
    plt.vlines(4+duration, 0.45, 0.6,linestyle ='--', color='r')
    plt.plot(t,y_8[:,0], label = r'$t=8$', color='green')
    plt.plot(t,y_8[:,2],alpha=alpha,linestyle ='--', color='green')

    plt.vlines(8, 0.15, 0.3, linestyle ='-', color='green')    
    plt.vlines(8+duration, 0.15, 0.3,linestyle ='--', color='green')

    plt.xlabel(r"$t$")
    plt.ylabel(r"incidence")
    plt.title(r"$S(t)$ (continuous) and $R(t)$ (dashed)")    
    '''