#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 10:48:44 2021

@author: fra
"""

import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def  contourplot(fig,axes,indd,xvector, yvector, xlabel, ylabel, toplot, labels, save,numberoflines,fmt,cmap, levels,lines='on'):
    '''contourplots'''
    X,Y = np.meshgrid(xvector,yvector)
    
    
    index = 0
    for Rinfty,ax in zip(toplot,axes):
        '''
        if index<2:
            #fig = plt.figure(figsize=(2,2))
            #plt.subplots_adjust(left=0.25, bottom=0.18, right=0.94, top=0.94, wspace=0, hspace=0)
        else:
            #fig = plt.figure(figsize=(2.4,2))
            #plt.subplots_adjust(left=0.25, bottom=0.18, right=0.9, top=0.94, wspace=0, hspace=0)

        ax = fig.add_subplot(111)
        '''
        #levels = np.arange(np.min(Rinfty)-0.2*np.min(Rinfty), np.max(Rinfty)+0.2*np.min(Rinfty), (np.max(Rinfty)-np.min(Rinfty))/numberoflines )
        print(np.min(Rinfty),np.max(Rinfty))
        Z =np.reshape(Rinfty,X.shape).T
        #norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=abs(Z).min())
        #cset1=ax.contourf(X, Y, Z, levels, norm=norm,
         #           cmap=cm.get_cmap(cmap, len(levels) - 1))
        #cset2 = ax.contour(X, Y, Z, cset1.levels, colors='k')
        #ax.clabel(cset2, levels[:-1],fmt=fmt,
        #  fontsize=10)
        c = ax.pcolor(X,Y, Z,cmap = cmap,vmin = np.min(toplot), vmax = np.max(toplot))
        ax.set_xlabel(xlabel,fontsize = 11, labelpad=-1)
        ax.set_ylabel(ylabel,fontsize = 11, labelpad=-1)
        ax.tick_params(axis='both', labelsize=10)        
        if indd ==0:
            left, bottom, width, height = ax.get_position().bounds
            cax = fig.add_axes([left+0.85*width, bottom, 0.1, height])
            cax.axis("off")
            fig.colorbar(c,ax=cax, pad=-30)      
       

        ax.text(0.9, 0.05, labels[index], transform=ax.transAxes, size = 11, color='k')
        ax.tick_params(axis='both', labelsize=12)
    
        ax.set_xlabel(xlabel,fontsize = 11, labelpad =-1)
        ax.set_ylabel(ylabel,fontsize = 11, labelpad=-3)

        #Z =np.reshape(Rinfty,X.shape).T
        #norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=abs(Z).min())

        #cset1=ax.contourf(X, Y, Z, levels, norm=norm,
        #             cmap=cm.get_cmap(cmap, len(levels) - 1))
        if lines != 'off':
            cset2 = ax.contour(X, Y, Z, levels, colors='k', linestyles='--')
            ax.clabel(cset2, levels,fmt=fmt,fontsize=10, colors='k')
        else:
            levels = np.arange(np.min(Rinfty)-0.2*np.min(Rinfty), np.max(Rinfty)+0.2*np.min(Rinfty), (np.max(Rinfty)-np.min(Rinfty))/numberoflines )
            cset2 = ax.contour(X, Y, Z, levels, colors='k', linestyles='--')
            
        #ax.set_zlabel(zlabel,fontsize = 8)
        #ax.view_init(rotateview[0], rotateview[1])
        #plt.savefig("figures/%s.eps"%save[index],format='eps')
        index +=1
        
if __name__=="__main__":
    '''
    It accepts 2 arguments, file_name which is the name of the file where
    all the results of metapop_strats are stored, and the var that 
    was fixed in that run (see metapop_strats).```````
    
    The it generates the contourplots accordingly (see figures in the paper)
    and saves them in figures/name

    the name is given by the following:
        
        fixed_var_9_string1_string2
    
    where fixed_var is either c or duration,
    string1 is either R_infty, I_max, t_avg, and string_2 is either 
    local/local_to_global/global.
    
    '''
    
    
    
    ''' file name is the seed of the random gen '''
    
    #file_name = sys.argv[1]
    file_name = "30_c"
    #R_0 matrix
    fixed_var = "c"
                  
    fig, ax = plt.subplots(1,2, figsize=(6,3))
    plt.subplots_adjust(left=0.08, bottom=0.15, right=0.8, top=0.9, wspace=0.3, hspace=0.3)


    ngroups = 9
    gamma = [1.0]*ngroups
    size = 100
    df = pd.read_csv('contourplots/%s.csv'%file_name)
    R_sub = np.array(df['R_sub'].tolist())
    I_sub = np.array(df['I_sub'].tolist())
    t_sub = np.array(df['t_sub'].tolist())
    R_subpop = np.array(df['R_subpop'].tolist())
    I_subpop = np.array(df['I_subpop'].tolist())
    t_subpop = np.array(df['t_subpop'].tolist())
    R_glob = np.array(df['R_glob'].tolist())
    I_glob = np.array(df['I_glob'].tolist())
    t_glob = np.array(df['t_glob'].tolist())
    '''
     strength of the intervention (i.e. during int, new R_0 = c*old R_0)
    '''
    c = [0.8]*ngroups
    tauf = 100
    
    
    interventionthreshold = np.linspace(0.025,1, size)
    interventionduration = np.linspace(0.5, 10, size)
    
    labels = (r'$\mathbf{A}$',r'$\mathbf{B}$', r'$\mathbf{C}$')
    xlabel = r"$Threshold$"
    ylabel = r"$Duration$"
    toplot = (I_glob,I_sub, I_subpop)
    levels=[0.07,0.08,0.09,0.11]
    fmt =	{
 
                  0.05: r"$0.05$",
                  0.07: r"$0.07$",
                  0.08: r"$0.08$",
                  0.09: r"$0.09$",
                  0.11: r"$0.11$"}
    
    save=['contourplots/%s_9_R_infty_local'%file_name, 'contourplots/%s_9_R_infty_localtoglobal'%file_name, 'contourplots/%s_9_R_infty_global'%file_name]
    axes = [ax[0]]
    
    contourplot(fig, axes,-1,interventionthreshold,interventionduration, xlabel, ylabel, toplot, labels, save, numberoflines=4, fmt = fmt,cmap = 'coolwarm',levels=levels,lines='on')

    file_name = "30_duration"
    #R_0 matrix
    fixed_var = "duration"
    df = pd.read_csv('contourplots/%s.csv'%file_name)
    R_sub = np.array(df['R_sub'].tolist())
    I_sub = np.array(df['I_sub'].tolist())
    t_sub = np.array(df['t_sub'].tolist())
    R_subpop = np.array(df['R_subpop'].tolist())
    I_subpop = np.array(df['I_subpop'].tolist())
    t_subpop = np.array(df['t_subpop'].tolist())
    R_glob = np.array(df['R_glob'].tolist())
    I_glob = np.array(df['I_glob'].tolist())
    t_glob = np.array(df['t_glob'].tolist())    
    interventionduration = [2]*ngroups
    cvector = np.linspace(0.1,0.9,size)
    interventionthreshold = np.linspace(0.025,1, size)
    toplot = (I_glob,I_sub, I_subpop)

    labels = (r'$\mathbf{B}$',r'$\mathbf{E}$', r'$\mathbf{F}$')
    xlabel = r"$Threshold$"
    ylabel = r"$c$"

    save=['contourplots/%s_9_R_infty_local'%file_name, 'contourplots/%s_9_R_infty_localtoglobal'%file_name, 'contourplots/%s_9_R_infty_global'%file_name]
    
    #toplot=(R_sub, R_subpop, R_glob)
    axes = [ax[1]]

    contourplot(fig, axes,0,interventionthreshold,cvector, xlabel, ylabel, toplot, labels, save,numberoflines=4, fmt = fmt,cmap='coolwarm',levels=levels)
  

    plt.savefig("fig8.tiff",dpi=600)
    plt.savefig("fig8.eps")
        
        
        