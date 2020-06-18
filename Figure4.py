#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 15:19:36 2020

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

from matplotlib import rc
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)



ngroups = 9
gamma = [1.0]*9
tauf = 35
betaij = np.loadtxt('mixing_baseline.txt', delimiter=',')

c =[0.5]*9

interventiontime = [1.1]*9
interventionduration = [4]*9

SIR = SIR_model(ngroups, gamma, tauf,betaij, betain=0, betaoff=0, seed=1)

y=SIR.sir_intervention( c, [1], interventiontime, interventionduration, nt = 3000, epsilon=0.01, intervention='subgroup_threshold')

#y[:ngroups] is the S_1(t)... S_n(t) susceptible populations evolution,
#y[ngroups:2*ngroups] "I(t)"
#y[2*ngroups:] "R(t)"

t = np.linspace(0,tauf,3000)
plt.close()
fig,ax = plt.subplots(3,3, figsize=(5.5,5.5), sharex=True,sharey = True)
ax = ax.ravel()
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.96, top=0.96, wspace=0.2, hspace=0.2)

#plot I(t)
for i,sub in enumerate(ax):
    #S(t)
    #sub.plot(t, y[:,i], color='b')
    #I(t)
    sub.plot(t, y[:,i+ngroups], color='r')
    #R(t)
    #sub.plot(t, y[:,i+2*ngroups], color='g')
    #intervention
    #sub.vlines(SIR.intervention_time[i], 0,np.max(y[:,i+ngroups]))
    sub.set_title("sub-population %d" %(i+1))

finalsize = np.sum(y[-1:,2*ngroups:])
ax[7].set_xlabel(r"$t$", size=11)
ax[3].set_ylabel(r"$I(t)$",size=11,labelpad=-2)

plt.savefig("9groupsepidemics_no_control.eps",format='eps', dpi=500)

#Figure 4