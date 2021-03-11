#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 10:03:58 2020

@author: ld288
"""

from Eulerclasssir import SIR_model
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from scipy.interpolate import interp1d
from scipy.integrate import quad
import sys
import pandas as pd

'''
The way this code works is analogous to figure 3, just extended to 9 groups
'''


def computetime(sol,tauf,nt, beta, ngroups):
    t = np.linspace(0,tauf,nt)
    dS = np.zeros((ngroups,nt-1))
    for i in range(ngroups):
        dS[i] = np.diff(sol[:,i])
        #print(sol[0,i])
    dS = np.sum(dS,axis=0)
    t_infection = np.sum(-t[:-1]*dS)  
    R_infty = np.sum(sol[-1,2*ngroups:])     
    return(t_infection/R_infty)







def figure0(gamma, beta, tauf, c, interventiontime, interventionduration, whichx ='c', whichy='time0', whichz='rinfty', typeofintervention='time', savename='fig1.eps',rotateview=(0,0)):
    '''
    Simply runs the euclassir code and returns R infty, avg time and I_peak
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
    
    t = np.linspace(0,tauf,250)

    if whichx == 'c':
        variable1 = c
        xlabel = r'strength of intervention'
    elif whichx =='time0':
        variable1 = interventiontime
        if typeofintervention =='time':
            xlabel = r'time of intervention'
        else:
            xlabel = r'$Threshold$'
    elif whichx == r'duration':
        variable1 = interventionduration
        xlabel = r'duration of intervention'
    else:
        print("what's on x??")
        return -1
    
    if whichy == 'c':
        variable2 = c
        ylabel = r'strength of intervention'
    elif whichy =='time0':
        variable2 = interventiontime
        if typeofintervention =='time':
            ylabel = r'time of intervention'
        else:
            ylabel = 'threshold of intervention'    
    elif whichy == r'duration':
        variable2 = interventionduration
        ylabel = r'$Duration$'
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
    int_time=[]

    for i,inputx in enumerate(variable1):
        for j,inputy in enumerate(variable2):
            SIR = SIR_model(ngroups, gamma, tauf,beta, betain=0, betaoff=0, seed=1)

            x =[inputx]*ngroups
            y=[inputy]*ngroups
            if whichx =='c' and whichy=='time0':
                sol=SIR.sir_intervention( x, 1, y, interventionduration, nt = 3000, epsilon=1e-2, intervention=typeofintervention)
            elif whichx =='c' and whichy=='duration':
                sol=SIR.sir_intervention( x, 1, interventiontime, y, nt = 3000, epsilon=1e-2, intervention=typeofintervention)
            elif whichx =='time0' and whichy=='c':
                sol=SIR.sir_intervention(y, 1, x, interventionduration, nt = 3000, epsilon=1e-2, intervention=typeofintervention)
            elif whichx =='time0' and whichy=='duration':
                sol=SIR.sir_intervention(c , 1, x, y, nt = 3000, epsilon=1e-2, intervention=typeofintervention)
            elif whichx =='duration' and whichy=='c':
                sol=SIR.sir_intervention( y, 1, interventiontime, x, nt = 3000, epsilon=1e-2, intervention=typeofintervention)
            elif whichx =='duration' and whichy=='time0':    
                sol=SIR.sir_intervention( c, 1, y, x, nt = 3000, epsilon=1e-2, intervention=typeofintervention)

            globalI = np.zeros_like(t)

            for k in range(len(t)):
                globalI[k] = np.sum(sol[k, 9:18])
            
            tmax.append(t[np.argmax(globalI)]/9.0)
            Imax.append(np.max(globalI)/9.0)
            Rinfty.append(np.sum(sol[-1,2*ngroups:])/9.0)
            int_time.append(computetime(sol,tauf,beta=beta, nt=3000, ngroups = 9))

    #x = np.repeat(variable1, len(variable2))
    #y = np.tile(variable2, len(variable1))
    return (Rinfty,Imax,int_time)


if __name__=="__main__":
 
    '''
     This function takes 2 arguments in. The first one, is the name (= the seed)
     of the file where to find the mixing matrix for the 9 metapopulation model
     
     The second one is a var which should be either 'c' or 'duration' and tells
     what should be fixed. If 'c', the intervention strength is fixed, while
     threshold and duration vary. If 'duration', duration is fixed and c and 
     intervention_strength vary.
     
     Then, it proceeds to implement on each combination of parameters the three
     strategies - subgropup_threshold, subgroup_to_pop_threshold, pop_threshold
     and it saves the file.
    '''
       
    
    
    
    ''' file name is the seed of the random gen '''
   
    
    #file_name = int(sys.argv[1])
    #fixed_var = str(sys.argv[2])
    #R_0 matrix
    #betaij = np.loadtxt('mat_realizations/%d.txt'%file_name)
    file_name ="mixing_matrix"
    fixed_var="c"
    betaij = np.loadtxt('mixing_baseline.txt',delimiter=',')
    ngroups = len(betaij)
    gamma = [1.0]*ngroups

    size = 50
    '''
     strength of the intervention (i.e. during int, new R_0 = c*old R_0)
    '''
    c = [0.8]*ngroups
    tauf = 50
        
    '''
        the next three runs are done with a fixed value of c and varying the threshold and the duration
    '''
    if fixed_var=='c':
        interventionthreshold = np.linspace(0.025,1, size)
        interventionduration = np.linspace(0.5, 10, size)
        
        # Simulations:
        (R_sub,I_sub,t_sub)= figure0(gamma,betaij,tauf,c,interventionthreshold,interventionduration,whichx ='time0', whichy='duration', whichz='rinfty', typeofintervention='subgroup_threshold', savename ='control_02', rotateview =(20,15))
        (R_subpop,I_subpop,t_subpop)= figure0(gamma,betaij,tauf,c,interventionthreshold,interventionduration,whichx ='time0', whichy='duration', whichz='rinfty', typeofintervention='subgroup_to_pop_threshold', savename ='control_02.eps', rotateview =(20,15))
        (R_glob,I_glob,t_glob)= figure0(gamma,betaij,tauf,c,interventionthreshold,interventionduration,whichx ='time0', whichy='duration', whichz='rinfty', typeofintervention='pop_threshold', savename ='control_02.eps', rotateview =(20,15))
        
       
        df = pd.DataFrame({'R_sub': R_sub,
                       'I_sub': I_sub, 
                       't_sub':t_sub,
                       'R_subpop': R_subpop,
                       'I_subpop': I_subpop, 
                       't_subpop':t_subpop,
                       'R_glob': R_glob,
                       'I_glob': I_glob, 
                       't_glob':t_glob,
                       })
    
        df.to_csv('contourplots/%s_c.csv'%file_name)

    elif fixed_var=='duration':  
        '''
            the next three runs are done with a fixed duration and varying the threshold and the strengths of intervention
        '''
        interventionduration = [2]*9
        cvector = np.linspace(0.1,0.9,size)
        interventionthreshold = np.linspace(0.025,1, size)
        
        # Simulations:
        (R_sub,I_sub,t_sub)= figure0(gamma,betaij,tauf,cvector,interventionthreshold,interventionduration,whichx ='time0', whichy='c', whichz='rinfty', typeofintervention='subgroup_threshold', savename ='control_02.eps', rotateview =(20,15))
        (R_subpop,I_subpop,t_subpop)= figure0(gamma,betaij,tauf,cvector,interventionthreshold,interventionduration,whichx ='time0', whichy='c', whichz='rinfty', typeofintervention='subgroup_to_pop_threshold', savename ='control_02.eps', rotateview =(20,15))
        (R_glob,I_glob,t_glob)= figure0(gamma,betaij,tauf,cvector,interventionthreshold,interventionduration,whichx ='time0', whichy='c', whichz='rinfty', typeofintervention='pop_threshold', savename ='control_02.eps', rotateview =(20,15))
        
        df = pd.DataFrame({'R_sub': R_sub,
                       'I_sub': I_sub, 
                       't_sub':t_sub,
                       'R_subpop': R_subpop,
                       'I_subpop': I_subpop, 
                       't_subpop':t_subpop,
                       'R_glob': R_glob,
                       'I_glob': I_glob, 
                       't_glob':t_glob,
                       })
    
        df.to_csv('contourplots/%s_duration.csv'%file_name)        

    else:
        print("Sry cannot fix this var")
    
    