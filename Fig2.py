#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 17:01:49 2020
@author: Francesco Di Lauro
@mail: F.Di-Lauro@sussex.ac.uk
Copyright 2020 Francesco Di Lauro. All Rights Reserved.
See LICENSE file for details
"""
from Eulerclasssir import SIR_model
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

from matplotlib import rc
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

ngroups=1
beta = [2.5]
gamma = [1]
c = [0.8]
tauf=30

#interventiontime = np.linspace(0.005,5, 1000)
interventionduration=[2]
t = np.linspace(0,tauf, 2000)
interventiontime = t

fig,ax = plt.subplots(3,2, figsize=(6,4.5))
ax = ax.ravel()
fig.subplots_adjust(wspace=0.42, hspace=0.2)
colors = ['r','b','g','pink','cyan', 'orange']


interventiontime = np.array([7, 4.8, 4.5, 3.6, 3, 1])

SIR = SIR_model(ngroups, gamma, tauf,beta, betain=0, betaoff=0, seed=1)

'''
We plot different scenarios where we vary the time of intervention,
hence intervention='time'
'''
for index,time in enumerate(interventiontime):
    SIR = SIR_model(ngroups, gamma, tauf,beta, betain=0, betaoff=0, seed=1)
    sol=SIR.sir_intervention( c, 0, [time], interventionduration, nt = len(t), epsilon=1e-3, intervention='time')   
    ax[1].plot(t, sol[:,0], color=colors[index], linestyle ='-')#,label =r'$t^\star = %d$'%int(time))
    #ax[1].vlines(time, 0, 1, color= colors[index])
    ax[3].plot(t, sol[:,1], color=colors[index], linestyle ='-')
    #ax[3].vlines(time, 0, np.max(sol[:,1]),color= colors[index])
    ax[5].plot(t, sol[:,1]+sol[:,2], color=colors[index], linestyle ='-')
    #ax[5].vlines(time, 0, np.max(sol[:,2]),color= colors[index])
    array_index = np.argwhere(t>time)[0][0]
    #print(array_index)
    Thresh = sol[array_index,1]+sol[array_index,2]
    print(interventiontime, Thresh)
    #print(Thresh)
    ax[0].axvline(Thresh, 0, 0.2, color= colors[index])
    ax[2].axvline(Thresh, 0, 0.2, color= colors[index])
    ax[4].axvline(Thresh, 0, 0.2, color= colors[index])
#ax[1].legend(loc='best')

sol=SIR.sir_intervention( c, 0, [2], interventionduration, nt = len(t), epsilon=1e-3, intervention='subgroup_threshold')  
ax[1].plot(t, sol[:,0], color='k', linestyle ='--',label =r'No intervention')

ax[3].plot(t, sol[:,1], color='k', linestyle ='--')
ax[5].plot(t, sol[:,1]+sol[:,2], color='k', linestyle ='--')

def computetime(sol,tauf,nt, beta, ngroups=1):
    t = np.linspace(0,tauf,nt)
    dS = np.zeros((ngroups,nt-1))
    for i in range(ngroups):
        dS[i] = np.diff(sol[:,i])
        #print(sol[0,i])
    dS = np.sum(dS,axis=0)
    t_infection = np.sum(-t[:-1]*dS)
  
    R_infty = sol[-1,2]
     
    return(t_infection/R_infty)   
#For panels 0,2,4 you need to plot as a function of threshold:
thresholdv = np.linspace(0,1,100)
rinftyv = np.zeros_like(thresholdv)
Ipeak = np.zeros_like(thresholdv)
tpeak = np.zeros_like(thresholdv)
interventionduration=[2]
for threshindex,threshold in enumerate(thresholdv):
    SIR =SIR_model(ngroups, gamma, tauf,beta, betain=0, betaoff=0, seed=1)
    sol=SIR.sir_intervention( c, 0, [threshold], interventionduration, nt = len(t), epsilon=1e-3, intervention='subgroup_threshold')
    rinftyv[threshindex] = sol[-1,2]
    Ipeak[threshindex] = np.max(sol[:,1])
    tpeak[threshindex] = computetime(sol,30,2000,beta)

ax[0].plot(thresholdv,rinftyv,color='k')
ax[2].plot(thresholdv,Ipeak,color='k')
ax[4].plot(thresholdv,tpeak,color='k')
index+=3


'''
This bit is wrong
for sub,plot in zip(ax[::2],toplot):
    sub.plot(sol[:,1]+sol[:,2], plot, color= 'g')
    index+=1
'''
for i in range(6):

    ax[i].tick_params(axis='both', labelsize=11)
  
label = [r'$R(\infty)$',r"$S(t)$", r'$I_{{max}}$',r"$I(t)$",r'$\overline{t}$',r"$I(t)+R(t)$",]
num = [r'$\mathbf{A}$',r'$\mathbf{B}$',r'$\mathbf{C}$',r'$\mathbf{D}$',r'$\mathbf{E}$',r'$\mathbf{F}$']
posx =[0.03,19]*3
posy =[0.76,0.85,0.24,0.25,7.3,0.1]
for i,sub in enumerate(ax):
    sub.set_ylabel(label[i],size= 12, labelpad=0.2)
    if i%2 ==0:

            
        sub.set_xlim(0,0.9)
        sub.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        sub.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
        sub.set_xticklabels([])

        if i ==4:
            sub.set_xlabel(r"$Threshold$", size= 12)
            sub.set_xticklabels([r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$'])
            sub.set_ylim(4,8)
            sub.set_yticks([4,6,8])
            sub.set_yticklabels([r'$4.0$',r'$6.0$',r'$8.0$'])
        if i==0:
            sub.set_ylim(0.75,0.90)
            sub.set_yticklabels([r'$0.75$',r'$0.80$',r'$0.85$',r'$0.90$'])
            
        if i==2:
            sub.set_ylim(0.1,0.26)
            sub.set_yticks([0.1,0.15,0.20,0.25])
            sub.set_yticklabels([r'$0.10$',r'$0.15$',r'$0.20$',r'$0.25$'])
    else:
        sub.set_xticks([0, 5, 10, 15, 20])
        sub.set_xticklabels([])
        sub.set_xlim(0,22)
        sub.set_ylim(0,.9)
        sub.set_ylim(0,1)

        if i ==5:
            sub.set_xlabel(r"$t$",size= 12)
            sub.set_xticklabels([r'$0$',r'$5$',r'$10$',r'$15$',r'$20$'])
        else:
            if i ==3:
                sub.set_ylim(0,0.30)

        
    #sub.text(posx[i], posy[i], num[i], size= 11)
    sub.text(0.03, 0.35, num[i],  transform=sub.transAxes, size= 11)
    sub.tick_params(axis='both', labelsize=12)

plt.show()
plt.savefig("fig2.tiff",dpi=600)
plt.savefig("fig2.eps")