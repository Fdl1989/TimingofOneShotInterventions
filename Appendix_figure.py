#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 15:19:36 2020

@author: Francesco Di Lauro
@mail: F.Di-Lauro@sussex.ac.uk
Copyright 2020 Francesco Di Lauro. All Rights Reserved.
See LICENSE file for details
"""

from Joel_contourplot_func import *
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from matplotlib import rc
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
'''
This code produces contourplots for different variables.
NOTE, this works using mixing_baseline.txt, which is now moved at the appendix.
Please refer to the files matrix_generation and metapop_strats
'''
def setup_plot(X, Y, Z, xlabel, ylabel, label):
    numberoflines = 4
    levels = np.arange(np.min(Z)-0.2*np.min(Z), np.max(Z)+0.2*np.min(Z), (np.max(Z)-np.min(Z))/numberoflines )
    fig = plt.figure(figsize=(2,2))
    plt.subplots_adjust(left=0.25, bottom=0.18, right=0.94, top=0.94, wspace=0, hspace=0)
    ax = fig.add_subplot(111)
    norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=abs(Z).min())
    cmap = cm.viridis
    
    cset1=ax.contourf(X, Y, Z, levels, norm=norm,
                        cmap=cm.get_cmap(cmap, len(levels) - 1))
    cset2 = ax.contour(X, Y, Z, cset1.levels, colors='k')
    ax.clabel(cset2, levels[:-1],fmt="%.2f",
            fontsize=10)
    ax.set_xlabel(xlabel,fontsize = 11, labelpad=-1)
    ax.set_ylabel(ylabel,fontsize = 11, labelpad=-1)
    ax.tick_params(axis='both', labelsize=10)        
    

    ax.text(0.05, 0.05, label, transform=plt.gca().transAxes, size = 10, color='r')
    ax.tick_params(axis='both', labelsize=12)
    
    ax.set_xlabel(xlabel,fontsize = 11, labelpad =-1)
    ax.set_ylabel(ylabel,fontsize = 11, labelpad=-3)
        #ax.set_zlabel(zlabel,fontsize = 8)
        #ax.view_init(rotateview[0], rotateview[1])
        

ngroups = 9
gamma = [1.0]*ngroups
c = [0.8]*9
tauf = 100
betaij = np.loadtxt('mixing_baseline.txt', delimiter=',')
'''
   Test of my code
'''
#interventiontime = [10]*9
interventionduration = [4]*9
#tot_I,peak_I,avg_time,S_control,I_control,R_control = epid_control(gamma,betaij,tauf,c,interventiontime,interventionduration,True)


int_times_strength = np.linspace(0,100,601)
strengths_strength = np.linspace(0,1,31)
X,Y_strength = np.meshgrid(int_times_strength, strengths_strength)
#print(X)
#print(Y)
global_size_strength = 0*X
max_size_strength = 0*X
tot_I_strength=0*X
peak_I_strength = 0*X
avg_time_strength = 0*X
for i, int_time in enumerate(int_times_strength):
    for j, strength in enumerate(strengths_strength):
        #print(int_time, strength)
        interventiontime = [int_time]*ngroups
        c = [strength]*ngroups
        #print(interventiontime, strength)
        tot_I_strength[j,i],peak_I_strength[j,i],avg_time_strength[j,i],S_control,I_control,R_control = epid_control(gamma,betaij,tauf,c,interventiontime,interventionduration,False)
        
        global_size_strength[j,i] = (sum(I_control)+sum(R_control))/ngroups
        max_size_strength[j,i] = max(1-S for S in S_control)


#print(tot_I)
#print(peak_I)
#plt.contourf()
#plt.show()
setup_plot(global_size_strength, Y_strength, tot_I_strength, r"$Threshold$", r"$c$", '$(f)$')
#plt.savefig('joel_contour
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.xticks([0,0.2,0.4,0.6, 0.8])
plt.savefig("figures/contourplots_9families/9_R_infty_global_strength.png", dpi=200)

setup_plot(max_size_strength, Y_strength, tot_I_strength,  r"$Threshold$", r"$c$", '$(e)$')
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.xticks([0,0.2,0.4,0.6, 0.8])
plt.savefig("figures/contourplots_9families/9_R_infty_local2global_strength.png", dpi=200)


plt.figure()
setup_plot(global_size_strength, Y_strength, avg_time_strength,  r"$Threshold$", r"$c$", '$(f)$')
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.xticks([0,0.2,0.4,0.6, 0.8])
plt.savefig("figures/contourplots_9families/9_ave_time_global_strength.png", dpi=200)

plt.figure()
setup_plot(max_size_strength, Y_strength, avg_time_strength,  r"$Threshold$", r"$c$", '$(e)$')
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.xticks([0,0.2,0.4,0.6, 0.8])
plt.savefig("figures/contourplots_9families/9_ave_time_local2global_strength.png", dpi=200)



plt.figure()
setup_plot(global_size_strength, Y_strength, peak_I_strength,  r"$Threshold$", r"$c$", '$(f)$')
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.xticks([0,0.2,0.4,0.6, 0.8])
plt.savefig("figures/contourplots_9families/9_peak_I_global_strength.png", dpi=200)

plt.figure()
setup_plot(max_size_strength, Y_strength, peak_I_strength,  r"$Threshold$", r"$c$", '$(e)$')
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.xticks([0,0.2,0.4,0.6, 0.8])
plt.savefig("figures/contourplots_9families/9_peak_I_local2global_strength.png", dpi=200)


ngroups = 9
gamma = [1.0]*ngroups
c = [0.8]*9
tauf = 100
betaij = np.loadtxt('mixing_baseline.txt', delimiter=',')
'''
   Test of my code
'''
#interventiontime = [10]*9
interventionduration = [4]*9
#tot_I,peak_I,avg_time,S_control,I_control,R_control = epid_control(gamma,betaij,tauf,c,interventiontime,interventionduration,True)


int_times = np.linspace(0,100,601)
durations = np.linspace(0,10,31)
X,Y_duration = np.meshgrid(int_times, durations)
#print(X)
#print(Y)
global_size_duration = 0*X
max_size_duration = 0*X
tot_I_duration=0*X
peak_I_duration = 0*X
avg_time_duration = 0*X
c = [0.8]*ngroups
for i, int_time in enumerate(int_times):
    for j, duration in enumerate(durations):
        #print(int_time, duration)
        interventionduration = [duration]*9
        interventiontime = [int_time]*ngroups
        
        #print(interventiontime, duration)
        tot_I_duration[j,i],peak_I_duration[j,i],avg_time_duration[j,i],S_control,I_control,R_control = epid_control(gamma,betaij,tauf,c,interventiontime,interventionduration,False)
        
        global_size_duration[j,i] = (sum(I_control)+sum(R_control))/ngroups
        max_size_duration[j,i] = max(1-S for S in S_control)


#print(tot_I)
#print(peak_I)
#plt.contourf()
#plt.show()
setup_plot(global_size_duration, Y_duration, tot_I_duration, r"$Threshold$", r"$Duration$", '$(c)$')
#plt.savefig('joel_contour
plt.yticks([0, 2.5, 5.0, 7.5, 10])
plt.xticks([0,0.2,0.4,0.6, 0.8])
plt.savefig("figures/contourplots_9families/9_R_infty_global_duration.png", dpi=200)

setup_plot(max_size_duration, Y_duration, tot_I_duration,  r"$Threshold$", r"$Duration$", '$(b)$')
plt.yticks([0, 2.5, 5.0, 7.5, 10])
plt.xticks([0,0.2,0.4,0.6, 0.8])
plt.savefig("figures/contourplots_9families/9_R_infty_local2global_duration.png", dpi=200)




plt.figure()
setup_plot(global_size_duration, Y_duration, avg_time_duration,  r"$Threshold$", r"$Duration$", '$(c)$')
plt.yticks([0, 2.5, 5.0, 7.5, 10])
plt.xticks([0,0.2,0.4,0.6, 0.8])
plt.savefig("figures/contourplots_9families/9_ave_time_global_duration.png", dpi=200)

plt.figure()
setup_plot(max_size_duration, Y_duration, avg_time_duration,  r"$Threshold$", r"$Duration$", '$(b)$')
plt.yticks([0, 2.5, 5.0, 7.5, 10])
plt.xticks([0,0.2,0.4,0.6, 0.8])
plt.savefig("figures/contourplots_9families/9_ave_time_local2global_duration.png", dpi=200)



plt.figure()
setup_plot(global_size_duration, Y_duration, peak_I_duration,  r"$Threshold$", r"$Duration$", '$(c)$')
plt.yticks([0, 2.5, 5.0, 7.5, 10])
plt.xticks([0,0.2,0.4,0.6, 0.8])
plt.savefig("figures/contourplots_9families/9_peak_I_global_duration.png", dpi=200)

plt.figure()
setup_plot(max_size_duration, Y_duration, peak_I_duration,  r"$Threshold$", r"$Duration$", '$(b)$')
plt.yticks([0, 2.5, 5.0, 7.5, 10])
plt.xticks([0,0.2,0.4,0.6, 0.8])
plt.savefig("figures/contourplots_9families/9_peak_I_local2global_duration.png", dpi=200)

