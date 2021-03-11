#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 09:31:01 2021

@author: fra
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 20:08:19 2020
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
from scipy.interpolate import interp1d
from scipy.integrate import quad

from matplotlib import rc
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)




'''
code for figure 1
'''
def computetime(sol,tauf,nt, beta):
    t = np.linspace(0,tauf,nt)
    dSdt = - beta[0]*sol[:,0]*sol[:,1]    
    dSdtfun = interp1d(t,dSdt)
    t_infection = quad(lambda x: -x*dSdtfun(x),0,t[-1])
    R_infty = sol[-1,-1]
    #print(t_infection[0]/R_infty)
    return t_infection[0]/R_infty
    





tauf = 15
ngroups = 1
beta =[2]*1
gamma = [1]

t = np.linspace(0,tauf,3000)
c =[0.5]

interventiontime = [1.1]
interventionduration = [4]

fig,ax = plt.subplots(1,2, figsize=(6,3.5), sharey = True)
ax = ax.ravel()
plt.subplots_adjust(wspace=0.1, hspace=0)

SIR = SIR_model(ngroups, gamma, tauf,beta, betain=0, betaoff=0, seed=1)

y=SIR.sir_intervention( c, [0], interventiontime, interventionduration, nt = 3000, epsilon=0.01, intervention='subgroup_threshold')

infectiontime=computetime(y,tauf,3000,beta)
ax[0].plot(t,y[:,0], color='b')
ax[0].plot(t,y[:,1], color='r')
ax[0].plot(t,y[:,2], color='g')
ax[0].hlines(np.max(y[:,1]),0,15, linestyle='--', color='k')
ax[0].vlines(infectiontime, 0, np.max(y[:,0]), linestyle='--', color='k')

ax[0].text(13, 0.5, r'$\mathbf{A}$', size = 12)

ax[0].set_ylabel(r"$S(t), I(t), R(t)$", size=10)
ax[0].set_xlabel(r"$t$", size=10)
ax[0].set_xlim(0,tauf)
ax[0].set_ylim(0,1)


#find attack rate:
import scipy.optimize as optimize

def func(x):
    return  1 - (1-0.01)*np.exp(-2*x)

R_0 = optimize.fixed_point(func,0.5)
ax[0].hlines(R_0,0,15, linestyle='--', color='g')
#ax[0].text(7,R_0+0.05, r"$R(\infty)$", size =11)



def func(x):
    return  1 - (1-0.01)*np.exp(-4*x)

R_0 = optimize.fixed_point(func,0.5)
ax[1].hlines(R_0,0,15, linestyle='--', color='g')
#ax[1].text(5,R_0+0.05, r"$R(\infty)$", size=11)




tauf = 10
t = np.linspace(0,tauf,3000)

beta =[4]*1
SIR = SIR_model(ngroups, gamma, tauf,beta, betain=0, betaoff=0, seed=1)

y=SIR.sir_intervention( c, [0], interventiontime, interventionduration, nt = 3000, epsilon=0.01, intervention='subgroup_threshold')
infectiontime=computetime(y,tauf,3000,beta)

ax[1].plot(t,y[:,0], color='b')
ax[1].plot(t,y[:,1], color='r')
ax[1].plot(t,y[:,2], color='g')
ax[1].set_xlabel(r"$t$", size=10)
ax[1].text(0.9, 0.5, r'$\mathbf{B}$', transform=plt.gca().transAxes, size = 12)
ax[1].hlines(np.max(y[:,1]),0,10, linestyle='--', color='k')
ax[1].vlines(infectiontime, 0, np.max(y[:,0]), linestyle='--', color='k')
ax[1].set_xlim(0,tauf)
ax[1].set_ylim(0,1)

plt.savefig("fig1.tiff",dpi=600)
plt.savefig("fig1.eps")