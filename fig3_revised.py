#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 10:03:58 2020

@author: Francesco Di Lauro
@mail: F.Di-Lauro@sussex.ac.uk
Copyright 2020 Francesco Di Lauro. All Rights Reserved.
See LICENSE file for details

"""

from Eulerclasssir import SIR_model
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from scipy.integrate import quad
from scipy.interpolate import interp1d
import pylab
from matplotlib import rc
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
'''
Might be useful to make a command that allows us to specify one of the variables, and a range for the other two and then generate the corresponding figure.
'''

def computetime(sol,tauf,nt, beta, ngroups=1):
    '''
    This function computes the average time of infection by summing
    - t dS/dt dt from 0 to infty
    '''
    t = np.linspace(0,tauf,nt)
    dS = np.zeros((ngroups,nt-1))
    for i in range(ngroups):
        dS[i] = np.diff(sol[:,i])
        #print(sol[0,i])
    dS = np.sum(dS,axis=0)
    t_infection = np.sum(-t[:-1]*dS)
  
    R_infty = np.sum(sol[-1,2*ngroups:])   
     
    return(t_infection/R_infty)    

def computetime_interpol(sol,tauf,nt,beta,ngroups=1):
    '''
    This function computes the average time of infection by interpolation, integrating
    - t dS/dt dt from 0 to infty
    '''
    t = np.linspace(0,tauf,nt)
    dS = np.zeros((ngroups,nt-1))
    dt = t[1]-t[0]
    for i in range(ngroups):
        dS[i] = np.diff(sol[:,i])
        #print(sol[0,i])     
    dS = np.sum(dS,axis=0)
    dSdt = dS/dt
    function = interp1d(t[:-1], dSdt)
    integrand =  lambda x: -x*function(x)
    integral = quad(integrand, 0, t[-2])[0]
    return integral/np.sum(sol[-1,2*ngroups:])


def simulation(gamma, beta, tauf, c, interventiontime, interventionduration, whichx ='c', whichy='time0', whichz='rinfty', typeofintervention='time'):
    '''
    gamma: recovery
    beta: infection
    tauf: finaltime
    c: intervention
    interventiontime: either time at which interventionstarts or threshold (depending on typeofintervention)
    interventionduration: duration of control
    whichx, whichy, whichz : what to plot on x,y,z
    typeofintervention: see interventiontime
    savename: name of the figure to save
    rotateview: angle of rotation of figure before saving
    
    '''
    
    t = np.linspace(0,tauf,2000)

    if whichx == 'c':
        variable1 = c
        xlabel = r'strength of intervention'
    elif whichx =='time0':
        variable1 = interventiontime
        if typeofintervention =='time':
            xlabel = r'$t^\star$'
        else:
            xlabel = r'$Threshold$'
    elif whichx == r'duration':
        variable1 = interventionduration
        xlabel = r'$Duration$'
    else:
        print("what's on x??")
        return -1
    
    if whichy == 'c':
        variable2 = c
        ylabel = r'$c$'
    elif whichy =='time0':
        variable2 = interventiontime
        if typeofintervention =='time':
            ylabel = r'time of intervention'
        else:
            ylabel = 'threshold of intervention'    
    elif whichy == r'duration':
        variable2 = interventionduration
        ylabel = r'$Duration$'
    elif whichy ==   'beta':
         variable2 = beta
         ylabel = r'$\mathcal{R}_0$'
       
    else:
        print("what's on y??")
        return -2
    
    if whichx == whichy:
        print("nasty 3d plot? x can't be equal to y")
        return -3
    
    #contour plot

    #fig = plt.figure(figsize=(4.5,4.5))
    #ax = fig.add_subplot(111, projection='3d')
    
    Imax = np.zeros(len(variable1)*len(variable2))
    tmax = np.zeros_like(Imax)
    Rinfty = np.zeros_like(Imax)
    interv_happ = np.zeros_like(Imax)
    S_hit = np.zeros_like(Imax)
    int_time = np.zeros_like(Imax)
    
    stupid_iterator=0
    for inputx in variable1:
        for inputy in variable2:
            SIR = SIR_model(ngroups, gamma, tauf,beta, betain=0, betaoff=0, seed=1)

            x =[inputx]*ngroups
            y=[inputy]*ngroups
            if whichx =='c' and whichy=='time0':
                sol=SIR.sir_intervention( x, 0, y, interventionduration, nt = 2000, epsilon=1e-3, intervention=typeofintervention)
            elif whichx =='c' and whichy=='duration':
                sol=SIR.sir_intervention( x, 0, interventiontime, y, nt = 2000, epsilon=1e-3, intervention=typeofintervention)
            elif whichx =='time0' and whichy=='c':
                sol=SIR.sir_intervention(y, 0, x, interventionduration, nt = 2000, epsilon=1e-3, intervention=typeofintervention)
            elif whichx =='time0' and whichy=='duration':
                sol=SIR.sir_intervention(c , 0, x, y, nt = 2000, epsilon=1e-3, intervention=typeofintervention)
            elif whichx =='duration' and whichy=='c':
                sol=SIR.sir_intervention( y, 0, interventiontime, x, nt = 2000, epsilon=1e-3, intervention=typeofintervention)
            elif whichx =='duration' and whichy=='time0':    
                sol=SIR.sir_intervention( c, 0, y, x, nt = 2000, epsilon=1e-3, intervention=typeofintervention)
            elif whichx =='time0' and whichy=='beta':
                #print(len(y), ngroups)
                SIR = SIR_model(ngroups, gamma, tauf,y, betain=0, betaoff=0, seed=1)

                sol=SIR.sir_intervention( c, 0, x, interventionduration, nt = 2000, epsilon=1e-3, intervention=typeofintervention)
            # the peak happens after the intervention is stopped, at the time of intervention, or before the intervention
            time = np.linspace(0,tauf,2000)
            indexpos = np.argmax(sol[:,1])
            tpeak = time[indexpos]
            if whichy=='beta':
                r0 = inputy/gamma[0]
                
            else:
                r0 = beta[0]/gamma[0]    
            try:
                time_end_inf = np.argwhere(time>SIR.intervention_time +SIR.duration)[0][0]
                

                #print(time[time_end_inf],tpeak)
                
                if sol[time_end_inf,0]*r0  < 1:
                    S_hit[stupid_iterator] = 0

                else:
                    S_hit[stupid_iterator] = 1
            except:
                if sol[-1,0] < 1/r0:
                    S_hit[stupid_iterator] = 0
                else:
                    S_hit[stupid_iterator] = 1
                
            if tpeak > SIR.intervention_time +SIR.duration:
                interv_happ[stupid_iterator] = 1
            elif tpeak > SIR.intervention_time:
                interv_happ[stupid_iterator] = 0
            else:
                interv_happ[stupid_iterator] = -1
            

            
            #print(tpeak, SIR.intervention_time)
            Imax[stupid_iterator] = np.max(sol[:,1])
            tmax[stupid_iterator] = SIR.intervention_time[0]
            Rinfty[stupid_iterator] = np.sum(sol[-1,ngroups:])
            int_time[stupid_iterator] = computetime(sol,tauf,beta=beta, nt=2000) 
            stupid_iterator+=1
    

    return (variable1,variable2,Rinfty,Imax,int_time, S_hit,interv_happ,xlabel,ylabel, whichz)




def figure(fig,ax,variable1,variable2,Rinfty,Imax,int_time,S_hit,interv_happ, xlabel, ylabel, whichz, savename='fig1.eps',label='j',colorbar=True):


    if whichz=='rinfty':
        '''
        if colorbar==True:
            fig = plt.figure(figsize=(2.4,2))
            plt.subplots_adjust(left=0.175, bottom=0.185, right=0.94, top=0.94, wspace=0, hspace=0)

        else:
            fig = plt.figure(figsize=(2,2))        
            plt.subplots_adjust(left=0.18, bottom=0.18, right=0.94, top=0.94, wspace=0, hspace=0)
        
        ax = fig.add_subplot(111)
        '''
        Z = np.reshape(Rinfty, (variable1.size,variable2.size))
        
        c = ax.pcolor(variable1,variable2, Z.T,cmap = 'plasma',vmin = 0, vmax =1)
        
        In =  np.reshape(S_hit, (variable1.size,variable2.size))
        X,Y = np.meshgrid(variable1,variable2)
        
        levels = np.array([0,1])
        cset2 = ax.contour(X, Y, In.T, levels, colors=['k'])
        fmt={}
        strs = [r'before','during','after']
        for l, s in zip(cset2.levels, strs):
            fmt[l] = s


        ax.set_ylim(np.min(variable2)+0.05, np.max(variable2))
        ax.set_xlabel(xlabel,fontsize = 11, labelpad=-1)
        ax.set_ylabel(ylabel,fontsize = 11, labelpad=-1)
        ax.tick_params(axis='both', labelsize=10)
       
        plt.rcParams["hatch.linewidth"] = 4



        #S(t^*+D) = 1/R_0 
        
        #levels = np.arange(np.min(Rinfty)-0.3, np.max(Rinfty)+0.3, (np.max(Rinfty)-np.min(Rinfty))/6  )
        #if label == 'g':
        #    levels = np.arange(np.min(Rinfty)-0.3, np.max(Rinfty)+0.3, (np.max(Rinfty)-np.min(Rinfty))/4  )
        levels = [0.3,0.5,0.7,0.8,0.9]    
        #Z =np.reshape(Rinfty,X.shape).T
        #norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=abs(Z).min())

        #cset1=ax.contourf(X, Y, Z, levels, norm=norm,
        #             cmap=cm.get_cmap(cmap, len(levels) - 1))
        cset2 = ax.contour(X, Y, Z.T, levels, colors='w', linestyles='--')
        ax.clabel(cset2, levels,fmt='%.1f',
          fontsize=10, colors='k')
        #ax.set_ylim(np.min(variable2)+0.05, np.max(variable2))
        #ax.text(0.86, 0.9, r'$(%s)$'%label, color='r', transform=plt.gca().transAxes, size = 12)

        #ax.tick_params(axis='both', labelsize=10)
       
        #plt.rcParams["hatch.linewidth"] = 4
        #rec1 = plt.Rectangle((0.75,0),1-0.75,6, facecolor="gray", alpha=0.6,
        #edgecolor="black", hatch=r"\\" )
        #ax.add_patch(rec1)
        
    elif whichz=='peaktime':
        '''
        if colorbar==True:
            fig = plt.figure(figsize=(2.4,2))
            plt.subplots_adjust(left=0.175, bottom=0.185, right=0.94, top=0.94, wspace=0, hspace=0)
        else:
            fig = plt.figure(figsize=(2,2))
            plt.subplots_adjust(left=0.18, bottom=0.18, right=0.94, top=0.94, wspace=0, hspace=0)
        
        ax = fig.add_subplot(111)
        '''
        Z = np.reshape(int_time, (variable1.size,variable2.size))
        
        c = ax.pcolor(variable1,variable2, Z.T,cmap = 'PuBu',vmin = 0, vmax = 16)
        
        
        
        ax.set_ylim(np.min(variable2)+0.05, np.max(variable2))
        ax.set_xlabel(xlabel,fontsize = 11, labelpad=-1)
        ax.set_ylabel(ylabel,fontsize = 11, labelpad=-2)
        ax.tick_params(axis='both', labelsize=10)

        plt.rcParams["hatch.linewidth"] = 4
 
        #rec1 = plt.Rectangle((0.75,0),1-0.75,6, facecolor="gray", alpha=0.6,
        #edgecolor="black", hatch=r"\\" )

        levels = [4,6,8,10,12,14]
        fmt =	{
              4: r"$4$",
              6: r"$6$",
              8: r"$8$",
              10: r"$10$",
              12: r"$12$",
              14: r"$14$"}
        #levels = np.arange(np.min(int_time)-1, np.max(int_time)+1, (np.max(int_time)-np.min(int_time))/5  )
        cset2 = ax.contour(variable1, variable2, Z.T, levels, colors='w', linestyles='--')

        if label=="I":
            ax.clabel(cset2, levels,fmt=fmt,
              fontsize=10, colors='k',rightside_up=True)     
        else:
            ax.clabel(cset2, levels,fmt=fmt,
              fontsize=10, colors='k',rightside_up=False)     


        '''
        levels = np.arange(np.min(int_time)-1, np.max(int_time)+1, (np.max(int_time)-np.min(int_time))/6  )

        cset2 = ax.contour(X, Y, Z.T, levels, colors='w', linestyles='--')
        ax.clabel(cset2, levels[:-1],fmt='%.2f',
          fontsize=10, colors='k')     
        
        fig = plt.figure(figsize=(2.5,2.5))
        plt.subplots_adjust(left=0.18, bottom=0.18, right=0.94, top=0.94, wspace=0, hspace=0)

        ax = fig.add_subplot(111, projection='3d')
        plt.subplots_adjust(left=0.12, bottom=0.15, right=0.94, top=0.94, wspace=0, hspace=0)
        ax.grid(False)

        Z =np.reshape(int_time,X.shape).T
        ax.plot_surface(X,Y,Z,cmap='viridis', edgecolor='none')
        #ax.set_ylim(np.min(variable2)+0.05, np.max(variable2))        
        #ax.text(0.86, 0.9, r'$(c)$', transform=plt.gca().transAxes, size = 12)
        ax.set_xlabel(r"$Threshold$",fontsize = 11, labelpad=-1)
        ax.set_ylabel(ylabel,fontsize = 11, labelpad=-1)
        ax.set_zlabel(r"$Mean$ $infection$ $time$",fontsize = 11, labelpad=-1)
        ax.tick_params(axis='both', labelsize=8)
        ax.text(0, 0, 8, r'$(%s)$'%label, color='red')
        ax.tick_params(axis='both', which='major', pad=-4)
        ax.view_init(30, 60)
        ax.set_zticks([5,10,15])   
        ax.set_zlim(0,max(int_time)+3)
        '''
    elif whichz == 'peakI':     
        '''
        if colorbar==True:
            fig = plt.figure(figsize=(2.4,2))
            plt.subplots_adjust(left=0.175, bottom=0.185, right=0.94, top=0.94, wspace=0, hspace=0)

        else:
            fig = plt.figure(figsize=(2,2))        

            plt.subplots_adjust(left=0.18, bottom=0.18, right=0.94, top=0.94, wspace=0, hspace=0)
        
        ax = fig.add_subplot(111)
        '''
        Z = np.reshape(Imax, (variable1.size,variable2.size))
        
        c = ax.pcolor(variable1,variable2, Z.T,cmap = 'viridis',vmin = 0, vmax=0.3)
        ax.set_ylim(np.min(variable2)+0.05, np.max(variable2))            
        In =  np.reshape(interv_happ, (variable1.size,variable2.size))
        X,Y = np.meshgrid(variable1,variable2)
        
        levels = np.array([-1,0,1])
        cset2 = ax.contour(X, Y, In.T, levels, colors=['r','y'])
        fmt={}
        strs = [r'before','during','after']
        for l, s in zip(cset2.levels, strs):
            fmt[l] = s
        #ax.clabel(cset2, levels,
        #  fontsize=10, fmt=[])        

       

        ax.tick_params(axis='both', labelsize=10)
        

        
        l#evels = np.arange(np.min(Imax)-0.2, np.max(Imax)+0.2, (np.max(Imax)-np.min(Imax))/5  )
        levels = [0.05,0.1,0.15,0.2,0.25,0.3]

        cset2 = ax.contour(X, Y, Z.T, levels, colors='w', linestyles='--')
        ax.clabel(cset2, levels,fmt='%.2f',
          fontsize=10, colors='k')  
        
    ax.set_xticks([0.2,0.4,0.6,0.8])  
    #ax.set_xticks([0.1,0.2,0.3,0.4])
    ax.set_xlim(0.05,1)
    ax.tick_params(axis='both', labelsize=10)
    if colorbar==True:
        left, bottom, width, height = ax.get_position().bounds
        cax = fig.add_axes([left+0.85*width, bottom, 0.1, height])
        cax.axis("off")
        fig.colorbar(c,ax=cax, pad=-30)     
        
    ax.text(0.86, 0.9, r'$\mathbf{%s}$'%label, color='k', size=12,transform=ax.transAxes)
    ax.set_xlabel(xlabel,fontsize = 11, labelpad=0)
    ax.set_ylabel(ylabel,fontsize = 11, labelpad=0)
    
    #ax.set_zlabel(zlabel,fontsize = 8)
    #ax.view_init(rotateview[0], rotateview[1])
    plt.savefig(savename,format='eps',pad_inches = 0)
    return ax
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.close()
ngroups = 1
beta =2.5
gamma = 1
tauf = 30
c = 0.8

t = np.linspace(0,tauf,2000)


beta = [beta]
gamma = [gamma]
c = [c]

'''
 First vary duration  with c and beta fixed
'''
fig, ax = plt.subplots(3,3, figsize=(6,6))
plt.subplots_adjust(left=0.08, bottom=0.08, right=0.85, top=0.9, wspace=0.25, hspace=0.25)

interventiontime = np.linspace(0.05,1, 100)
interventionduration = np.linspace(0.1, 6, 100)

#Data=simulation(gamma,beta,tauf,c,interventiontime,interventionduration,whichx ='time0', whichy='duration', whichz='rinfty', typeofintervention='subgroup_threshold')
#Data = np.array(Data,dtype=object)
savename ='SIR_onegroup_rinfty_duration_cmap.eps'
label='A'

#np.save("Data/fig3_t_d_r.npy", Data)
Data=np.load("Data/fig3_t_d_r.npy", allow_pickle=True)
figure(fig, ax[0,0],*Data,savename,label,colorbar=False)


#Data=simulation(gamma,beta,tauf,c,interventiontime,interventionduration,whichx ='time0', whichy='duration', whichz = 'peakI',typeofintervention='subgroup_threshold')
#Data = np.array(Data,dtype=object)
savename ='SIR_onegroup_peakI_duration_cmap.eps'
label='D'

#np.save("Data/fig3_t_d_p.npy", Data)
Data=np.load("Data/fig3_t_d_p.npy", allow_pickle=True)
figure(fig, ax[1,0],*Data,savename,label,colorbar=False)



#Data=simulation(gamma,beta,tauf,c,interventiontime,interventionduration,whichx ='time0', whichy='duration', whichz = 'peaktime',typeofintervention='subgroup_threshold')
#Data = np.array(Data,dtype=object)
savename ='SIR_onegroup_time_duration_cmap.eps'
label='G'

#np.save("Data/fig3_t_d_t.npy", Data)
Data=np.load("Data/fig3_t_d_t.npy", allow_pickle=True)
figure(fig, ax[2,0],*Data,savename,label,colorbar=False)

 

#Then vary $c$ holding $\beta$ and duration fixed in the next row. 

interventionduration = [4]
interventionc = np.linspace(0.2,0.9, 100)


#Data=simulation(gamma,beta,tauf,interventionc,interventiontime,interventionduration,whichx ='time0', whichy='c', whichz='rinfty', typeofintervention='subgroup_threshold')
#Data = np.array(Data,dtype=object)
savename ='SIR_onegroup_rinfty_c_cmap.eps'
label='B'

#np.save("Data/fig3_t_c_r.npy", Data)
Data=np.load("Data/fig3_t_c_r.npy", allow_pickle=True)
b=figure(fig, ax[0,1],*Data,savename,label,colorbar=False)


#Data=simulation(gamma,beta,tauf,interventionc,interventiontime,interventionduration,whichx ='time0', whichy='c', whichz = 'peakI',typeofintervention='subgroup_threshold')
#Data = np.array(Data,dtype=object)
savename ='SIR_onegroup_peakI_c_cmap.eps'
label='E'

#np.save("Data/fig3_t_c_p.npy", Data)
Data=np.load("Data/fig3_t_c_p.npy", allow_pickle=True)
e=figure(fig, ax[1,1],*Data,savename,label,colorbar=False)

#Data=simulation(gamma,beta,tauf,interventionc,interventiontime,interventionduration,whichx ='time0', whichy='c', whichz = 'peaktime',typeofintervention='subgroup_threshold')
#Data = np.array(Data,dtype=object)
savename ='SIR_onegroup_time_c_cmap.eps'
label='H'

#np.save("Data/fig3_t_c_t.npy", Data)
Data=np.load("Data/fig3_t_c_t.npy", allow_pickle=True)
h=figure(fig, ax[2,1],*Data,savename,label,colorbar=False)


#Then vary $\Ro$ from $1 to 4$, holding $c$ and $D$ in the final row 

interventionc = [0.8]
interventionduration = [4]
interventiontime = np.linspace(0.05,1, 100)
betav = np.linspace(1.05,4, 100)


#Data=simulation(gamma,betav,tauf,interventionc,interventiontime,interventionduration,whichx ='time0', whichy='beta', whichz='rinfty', typeofintervention='subgroup_threshold')
#Data = np.array(Data,dtype=object)
savename ='SIR_onegroup_rinfty_beta_cmap.eps'
label='C'

#np.save("Data/fig3_t_b_r.npy", Data)
Data=np.load("Data/fig3_t_b_r.npy", allow_pickle=True)
c=figure(fig, ax[0,2],*Data,savename,label,colorbar=True)


#Data=simulation(gamma,betav,tauf,interventionc,interventiontime,interventionduration,whichx ='time0', whichy='beta', whichz = 'peakI',typeofintervention='subgroup_threshold')
#Data = np.array(Data,dtype=object)
savename ='SIR_onegroup_peakI_beta_cmap.eps'
label='F'

#np.save("Data/fig3_t_b_p.npy", Data)
Data=np.load("Data/fig3_t_b_p.npy", allow_pickle=True)
figure(fig, ax[1,2],*Data,savename,label,colorbar=True)


#Data=simulation(gamma,betav,tauf,interventionc,interventiontime,interventionduration,whichx ='time0', whichy='beta', whichz = 'peaktime',typeofintervention='subgroup_threshold')
#Data = np.array(Data,dtype=object)
savename ='SIR_onegroup_time_beta_cmap.eps'
label='I'

#np.save("Data/fig3_t_b_t.npy", Data)
Data=np.load("Data/fig3_t_b_t.npy", allow_pickle=True)
i=figure(fig, ax[2,2],*Data,savename,label,colorbar=True)

ax[0,0].set_title(r"$Variable$ $Duration$")

ax[0,1].set_title(r"$Variable$ $c$")


ax[0,2].set_title(r"$Variable$ $R_0$")

plt.savefig("fig3.tiff",dpi=600)
plt.savefig("fig3.eps")




















#tauf=25
#t = np.linspace(0,tauf,1000)

'''
fig = plt.figure(figsize=(3.5,3.5))


SIR = SIR_model(ngroups, gamma, tauf,beta, betain=0, betaoff=0, seed=1)

threshold =[4]
duration=[23]
y=SIR.sir_intervention( c, 0, threshold, duration, nt = 1000, epsilon=1e-3, intervention='subgroup_threshold')   
plt.plot(t,y[:,2], color='k', linestyle='--')



interventiontime = [1.5, 3, 4.5]
interventionduraton = [2,4]
SIR = SIR_model(ngroups, gamma, tauf,beta, betain=0, betaoff=0, seed=1)

for time in interventiontime:
    for duration in interventionduraton:
        initial = [time]
        dur = [duration]
        y=SIR.sir_intervention( c, 0, initial, dur, nt = 1000, epsilon=1e-2, intervention='time')   
        if time ==1.5:
            plt.plot(t,y[:,2],color='b')
            
        elif time == 3:
            plt.plot(t,y[:,2],color='g')
        elif time == 4.5:
            plt.plot(t,y[:,2],color='r')
            
plt.vlines(1.5,0,0.85, colors='b', linestyle ='--')            
plt.vlines(3,0,0.85, colors='g', linestyle ='--')
plt.vlines(4.5,0,0.85, colors='r', linestyle ='--')        
        #plt.plot(t,y[:,1])
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
plt.xlim(0,tauf)
plt.ylim(0,1)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.xlabel(r"$t$")
plt.ylabel(r"$R(t)$", labelpad=-2)
plt.text(0.88, 0.92, r'$(a)$', transform=plt.gca().transAxes, size = 14)
plt.savefig("SIR_onegroup_controlexamplesR.eps", format='eps',dpi=400)
'''







'''


fig = plt.figure(figsize=(3.5,3.5))


SIR = SIR_model(ngroups, gamma, tauf,beta, betain=0, betaoff=0, seed=1)

threshold =[4]
duration=[23]
y=SIR.sir_intervention( c, 0, threshold, duration, nt = 1000, epsilon=1e-3, intervention='subgroup_threshold')   
plt.plot(t,y[:,1], color='k', linestyle='--')



interventiontime = [1.5, 3, 4.5]
interventionduraton = [2,4]
t = np.linspace(0,tauf,1000)
SIR = SIR_model(ngroups, gamma, tauf,beta, betain=0, betaoff=0, seed=1)

for time in interventiontime:
    for duration in interventionduraton:
        initial = [time]
        dur = [duration]
        y=SIR.sir_intervention( c, 0, initial, dur, nt = 1000, epsilon=1e-2, intervention='time')   
        if time ==1.5:
            plt.plot(t,y[:,1],color='b')
            
        elif time == 3:
            plt.plot(t,y[:,1],color='g')
        elif time == 4.5:
            plt.plot(t,y[:,1],color='r')
            
plt.vlines(1.5,0,0.25, colors='b', linestyle ='--')            
plt.vlines(3,0,0.25, colors='g', linestyle ='--')
plt.vlines(4.5,0,0.25, colors='r', linestyle ='--')        
        #plt.plot(t,y[:,1])
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.xlabel(r"$t$")
plt.ylabel(r"$I(t)$", labelpad = -2)
plt.xlim(0,tauf)
plt.ylim(0,0.25)
plt.text(0.88, 0.92, r'$(b)$', transform=plt.gca().transAxes, size = 14)
plt.savefig("SIR_onegroup_controlexamplesI.eps", format='eps',dpi=400)

#c = [0.5]
#interventionduration=np.linspace(0.5,4, 30)
#thresholds = np.linspace(0.01, 0.2, 30)
#figure0(gamma,beta,tauf,c,thresholds,interventionduration,whichx='time0', whichy='duration',  whichz = 'rinfty' , typeofintervention='threshold')
'''
'''

ngroups = 9
gamma = [1.0]*9
tauf = 200
betaij = np.loadtxt('mixing_baseline.txt', delimiter=',')

c =[0.5]*9

interventiontime = np.linspace(0.025,0.2, 30)
interventionduration = np.linspace(0.5, 4, 30)
figure0(gamma,betaij,tauf,c,interventiontime,interventionduration,whichx ='time0', whichy='duration', whichz='rinfty', typeofintervention='threshold', savename ='rinfty.eps', rotateview =(20,15))

x,y,z = np.loadtxt("scatter_data.txt", delimiter=',', unpack=True)

fig = plt.figure(figsize=(4.5,4.5))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x,y,z, c='r')
''' 