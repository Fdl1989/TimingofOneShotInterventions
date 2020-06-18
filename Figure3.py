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


def figure0(gamma, beta, tauf, c, interventiontime, interventionduration, whichx ='c', whichy='time0', whichz='rinfty', typeofintervention='time', savename='fig1.eps',label='j' ):
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
    
    Imax = []
    tmax = []
    Rinfty = []
    interv_happ=[]
    int_time=[]
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
            
            interv_happ.append(SIR.int_happ)
            Imax.append(np.max(sol[:,1]))
            tmax.append(SIR.intervention_time[0])
            Rinfty.append(np.sum(sol[-1,ngroups:]))
            int_time.append(computetime(sol,tauf,beta=beta, nt=2000)) 
            

    #x = np.repeat(variable1, len(variable2))
    #y = np.tile(variable2, len(variable1))
    X,Y = np.meshgrid(variable1,variable2)
    cmap = cm.viridis

    if whichz=='rinfty':
        fig = plt.figure(figsize=(2,2))
        plt.subplots_adjust(left=0.18, bottom=0.18, right=0.94, top=0.94, wspace=0, hspace=0)
        
        ax = fig.add_subplot(111)

        levels = np.arange(np.min(Rinfty)-0.3, np.max(Rinfty)+0.3, (np.max(Rinfty)-np.min(Rinfty))/6  )
        if label == 'g':
            levels = np.arange(np.min(Rinfty)-0.3, np.max(Rinfty)+0.3, (np.max(Rinfty)-np.min(Rinfty))/4  )
            
        Z =np.reshape(Rinfty,X.shape).T
        norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=abs(Z).min())

        cset1=ax.contourf(X, Y, Z, levels, norm=norm,
                     cmap=cm.get_cmap(cmap, len(levels) - 1))
        cset2 = ax.contour(X, Y, Z, cset1.levels, colors='k')
        ax.clabel(cset2, levels[:-1],fmt='%.2f',
          fontsize=10)
        ax.set_ylim(np.min(variable2)+0.05, np.max(variable2))
        ax.text(0.86, 0.9, r'$(%s)$'%label, color='r', transform=plt.gca().transAxes, size = 12)
        ax.set_xlabel(xlabel,fontsize = 11, labelpad=-1)
        ax.set_ylabel(ylabel,fontsize = 11, labelpad=-1)
        ax.tick_params(axis='both', labelsize=10)
       
        plt.rcParams["hatch.linewidth"] = 4
        rec1 = plt.Rectangle((0.75,0),1-0.75,6, facecolor="gray", alpha=0.6,
        edgecolor="black", hatch=r"\\" )
        #ax.add_patch(rec1)
        
    elif whichz=='peaktime':
        fig = plt.figure(figsize=(2,2))
        plt.subplots_adjust(left=0.19, bottom=0.18, right=0.95, top=0.94, wspace=0, hspace=0)

        ax = fig.add_subplot(111)
        levels = np.arange(np.min(int_time)-1, np.max(int_time)+1, (np.max(int_time)-np.min(int_time))/7 )
        if label == 'f':
            levels = np.arange(np.min(int_time)-1, np.max(int_time)+1, (np.max(int_time)-np.min(int_time))/4 )
        if label == 'i':
            levels = np.arange(np.min(int_time)-1, np.max(int_time)+1, (np.max(int_time)-np.min(int_time))/3 )
        
        Z =np.reshape(int_time,X.shape).T
        norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=abs(Z).min())

        cset1=ax.contourf(X, Y, Z, levels, norm=norm,
                     cmap=cm.get_cmap(cmap, len(levels) - 1))
        cset2 = ax.contour(X, Y, Z, cset1.levels, colors='k')
        ax.clabel(cset2, levels[:-1],fmt='%.1f',
          fontsize=10)
        ax.set_ylim(np.min(variable2)+0.05, np.max(variable2))
        ax.text(0.86, 0.9, r'$(%s)$'%label, color='r', transform=plt.gca().transAxes, size = 12)
        ax.set_xlabel(xlabel,fontsize = 11, labelpad=-1)
        ax.set_ylabel(ylabel,fontsize = 11, labelpad=-2)
        ax.tick_params(axis='both', labelsize=10)

        plt.rcParams["hatch.linewidth"] = 4
        #rec1 = plt.Rectangle((0.75,0),1-0.75,6, facecolor="gray", alpha=0.6,
        #edgecolor="black", hatch=r"\\" )
        '''
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
        
        fig = plt.figure(figsize=(2,2))
        plt.subplots_adjust(left=0.21, bottom=0.18, right=0.97, top=0.94, wspace=0, hspace=0)

        ax = fig.add_subplot(111)

        levels = np.arange(np.min(Imax)-0.1, np.max(Imax)+0.05, (np.max(Imax)-np.min(Imax))/3  )

        Z =np.reshape(Imax,X.shape).T
        norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=abs(Z).min())
        cset1=ax.contourf(X, Y, Z, levels, norm=norm,
                     cmap=cm.get_cmap(cmap, len(levels) - 1))
        cset2 = ax.contour(X, Y, Z, cset1.levels, colors='k')
        ax.clabel(cset2,fmt='%.2f',
          fontsize=10)
        ax.set_ylim(np.min(variable2)+0.05, np.max(variable2))            
        ax.text(0.86, 0.9,r'$(%s)$'%label,color='r', transform=plt.gca().transAxes, size = 12)
       
        ax.set_xlabel(xlabel,fontsize = 11, labelpad=-1)
        ax.set_ylabel(ylabel,fontsize = 11, labelpad=-1)
    
        ax.tick_params(axis='both', labelsize=10)
        
    ax.set_xticks([0.2,0.4,0.6,0.8])  
    #ax.set_xticks([0.1,0.2,0.3,0.4])
    ax.set_xlim(0.05,1)
    ax.tick_params(axis='both', labelsize=10)
    #ax.set_zlabel(zlabel,fontsize = 8)
    #ax.view_init(rotateview[0], rotateview[1])
    plt.savefig(savename,format='eps',pad_inches = 0, dpi=400)
    

    return (Rinfty,Imax,int_time)




'''
SIR = SIR_model(ngroups, gamma, tauf,beta, betain=0, betaoff=0, seed=1)

threshold =[4]
duration=[23]
y=SIR.sir_intervention( c, 0, threshold, duration, nt = 1000, epsilon=1e-3, intervention='subgroup_threshold')   
fig = plt.figure(figsize=(3.5,3.5))
plt.plot(t,y[:,0], color= 'b',label=r'S')
plt.plot(t,y[:,1], color='r', label=r'I')
plt.plot(t,y[:,2], color='g',label=r'R')
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.xlabel(r'$t$', size=11)
plt.ylabel(r'$I(t)$', size=11)
plt.text(0.88, 0.9, r'$(a)$', transform=plt.gca().transAxes, size = 11)
plt.legend(loc='best', fontsize=8)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
plt.savefig("SIR_onegroup_nocontrol.eps", format ='eps', dpi=400)
'''
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

interventiontime = np.linspace(0.05,1, 20)
interventionduration = np.linspace(0.1, 6, 20)


#figure0(gamma,beta,tauf,c,interventiontime,interventionduration,whichx ='time0', whichy='duration', whichz='rinfty', typeofintervention='subgroup_threshold',savename ='SIR_onegroup_rinfty_duration.eps',label='a')
#figure0(gamma,beta,tauf,c,interventiontime,interventionduration,whichx ='time0', whichy='duration', whichz = 'peakI',typeofintervention='subgroup_threshold', savename ='SIR_onegroup_peakI_duration.eps',label='b' )
#figure0(gamma,beta,tauf,c,interventiontime,interventionduration,whichx ='time0', whichy='duration', whichz = 'peaktime',typeofintervention='subgroup_threshold', savename ='SIR_onegroup_time_duration.eps',label='c' )
 
'''
Then vary $c$ holding $\beta$ and duration fixed in the next row. 
'''
interventionduration = [4]
interventionc = np.linspace(0.2,0.9, 20)
#figure0(gamma,beta,tauf,interventionc,interventiontime,interventionduration,whichx ='time0', whichy='c', whichz='rinfty', typeofintervention='subgroup_threshold',savename ='SIR_onegroup_rinfty_c.eps',label='d')
#figure0(gamma,beta,tauf,interventionc,interventiontime,interventionduration,whichx ='time0', whichy='c', whichz = 'peakI',typeofintervention='subgroup_threshold', savename ='SIR_onegroup_peakI_c.eps',label='e')
#figure0(gamma,beta,tauf,interventionc,interventiontime,interventionduration,whichx ='time0', whichy='c', whichz = 'peaktime',typeofintervention='subgroup_threshold', savename ='SIR_onegroup_time_c.eps',label='f' )
 
'''
Then vary $\Ro$ from $1 to 4$, holding $c$ and $D$ in the final row 
'''
interventionc = [0.8]
interventionduration = [4]
betav = np.linspace(1.2,4, 20)
figure0(gamma,betav,tauf,interventionc,interventiontime,interventionduration,whichx ='time0', whichy='beta', whichz='rinfty', typeofintervention='subgroup_threshold',savename ='SIR_onegroup_rinfty_beta.eps',label='g')
figure0(gamma,betav,tauf,interventionc,interventiontime,interventionduration,whichx ='time0', whichy='beta', whichz = 'peakI',typeofintervention='subgroup_threshold', savename ='SIR_onegroup_peakI_beta.eps',label='h' )
figure0(gamma,betav,tauf,interventionc,interventiontime,interventionduration,whichx ='time0', whichy='beta', whichz = 'peaktime',typeofintervention='subgroup_threshold', savename ='SIR_onegroup_time_beta.eps',label='i' )
 

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