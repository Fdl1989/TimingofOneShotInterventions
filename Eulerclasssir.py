
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 18:30:48 2020

@author: ld288
"""

from scipy.integrate import  odeint
import numpy as np
##############
from matplotlib import pyplot as plt

'''
The purpose of this code is to define a class (SIR_model) to do simulations of metapopulation SIR epidemics.

If you run this code you will get some tests

'''

class SIR_model():
    def __init__(self, ngroups, gamma, tauf, betamatrix, betain=0, betaoff=0, seed=1):
        '''
            Class to simulate a metapopulation SIR compartmental model on n groups with the following characteristics:
                1) Each group has contacts with the whole population, with rates of infections higher if the contact is in his group and lower otherwise
                2) Recovery rate is fixed for all of the groups
                3) Internal infection rates are distributed normally around a value B_0
                4) Inter-groups infection rates are distributed uniformly in [0,1/10 (max B_ii)]
                5) optionally, at a time T_0 we can lower the infection rate in a group by a percentage c, for a duration D
            
            Parameters:  
                ngroups: int, number of groups (default should be 9)
                beta0: float, main diagonal value of the infection matrix
                betamatrix: as input, otherwise random matrix will be generated (ignored if 0())
                betao: float, reduction of beta outside diagonal (ignored if 0)
                gamma[i]: float array, recovery rate (default should be 1)
                tauf: float, duration of the epidemic outbreak
                seed: float, seed of rng
        
        
        '''
        np.random.seed(seed)

        self.ngroups = ngroups
        self.maxir =0 
        self.gamma = gamma
        self.tauf = tauf
                
        if betain ==0 and betaoff ==0:
            self.betaij = betamatrix
        else:
            self.betain = betain
            self.betaoff = betaoff            
            self.prepare_beta()
       # if len(self.betaij) != self.ngroups:
       #     print("WARNING, len(betaij) different than ngroups")
    def prepare_beta(self):
        '''
                prepare the beta matrix in the following way. On the main diagonal you have a number close to betain, outside a number close to betaoff*betain.               
        '''
        self.betaij = np.zeros((self.ngroups,self.ngroups))
        for i in range(0,self.ngroups):
            self.betaij[i][i] = np.random.uniform(self.betain-0.5, self.betain+0.5)
            for j in range(0,self.ngroups):
                if j!=i:
                    self.betaij[i][j] = np.random.uniform(0,self.betaoff*(self.betain+0.5))
    
    
 
    def sir_intervention(self, c,  ingroup, threshold, duration, nt = 1000, epsilon =0.01, intervention='subgroup_threshold'):
        '''
            ODE for Joel system in groups, each group suppresses both beta in and beta out by a factor of (1-c) as soon as the number of infected reach a threshold
            
            Params:
                c: array(ngroups) reduction in percentage when intervetion happens, for each group
                ingroup: initial group to be infected
                threshold:array(ngroups) threshold at which intervention kicks in, for each group:
                    Note, if type of intervention is 'time' threshold is the :ime at which intervention happens.
                    If type of intervention is any of the others, such as subrgroup_threshold, then threshold refers to the
                    threshold when intervention happens.
                duration:array(ngroups) duration of intervention, for each group
                nt: timegrid
                epsilon: fraction of infected at time 0
                intervention: either threshold if you want the intervention to happen after a certain threshold, or time if you want them to happen at a given time independently
            Outputs:
                y, ndarray: solution to an odeint such that y[,0..ngroups] is the evolution of the susceptible in the groups, y[,ngrousp:2*ngroups] is the evolution of the infected and y[,2ngroups...3ngroups] is the evolution of recovered
        '''
        t = np.linspace(0,self.tauf,nt)
        
        self.c_red = c
        self.duration = duration
        self.threshold = threshold
        self.maxir=0
        self.intervention_time = np.array(self.ngroups*[self.tauf], dtype=float)
        self.intervention_index = np.array(self.ngroups*[self.tauf], dtype=int)
        self.int_happ = np.zeros(self.ngroups, dtype=bool)
        #set up initial condition, 0.01 percent of the population in group one is infected
        Sin = np.ones(self.ngroups)
        Sin[ingroup] = 1-epsilon
        Iin = np.zeros(self.ngroups)  
        Iin[ingroup] = epsilon
        Rin = np.zeros(self.ngroups)
        self.incond = np.hstack((Sin,Iin,Rin))       
        
        #call sir_interventionode to solve the ode
        if intervention =='subgroup_threshold':
            if threshold[0]>1:
                print("Warning: you choose a threshold and put a number >1. Maybe you want to use 'time' instead?")
            self.check = self.check_subgroupthreshold
        elif intervention =='time':
            self.check = self.check_time
        elif intervention =='pop_threshold':
            if threshold[0]>1:
                print("Warning: you choose a threshold and put a number >1. Maybe you want to use 'time' instead")
            self.check = self.check_popthreshold
        elif intervention =='subgroup_to_pop_threshold':
            if threshold[0]>1:
                print("Warning: you choose a threshold and put a number >1. Maybe you want to use 'time' instead")
            self.check = self.subgroup_to_pop_threshold            
        else:
            1/0
        y= self.sir_interventionode(self.incond, t)
        y = y[:-1]
        return y
    def beta_value(self,t):
        '''
           routine that updates the betamatrix according to the rates (modified in case of intervention)
        '''
        betanew = np.copy(self.betaij)
        if self.ngroups>1:
            for i in range(self.ngroups):
                #check that intervetion has kicked off and that is not already finished
                if self.int_happ[i] !=False and t<self.intervention_time[i] + self.duration[i]:
                    for j in range(self.ngroups):
                        betanew[i,j] = (1-self.c_red[i])*betanew[i,j]
                        #we don't want the diagonal to be reduced by a factor (1-c)^2
                        if j!=i:
                            betanew[j,i] = (1-self.c_red[i])*betanew[j,i]
        else:
            if self.int_happ[0] !=False and t<self.intervention_time[0] + self.duration[0]:
                betanew = (1-self.c_red[0])*betanew
        return betanew
    def check_subgroupthreshold(self, v,t):
        '''
         routine that keeps track of whether and when intervention has happened, each subgroup has its own threshold
        '''
        for i in range(self.ngroups):
            if v[i+self.ngroups]+v[i+2*self.ngroups] > self.threshold[i] and self.int_happ[i] == False:
                self.int_happ[i] = True
                self.intervention_time[i] = t

    def check_popthreshold(self, v,t):
        '''
            intervention kicks in as soon as the total cumulative number of infected reaches threshold
        '''
        totalthreshold=0
        for i in range(self.ngroups):
            totalthreshold +=v[i+self.ngroups]+v[i+2*self.ngroups] 
        if totalthreshold/self.ngroups > self.threshold[0] and self.int_happ[0] == False:
            maxir = 0
            
            for i in range(self.ngroups):
                self.int_happ[i] = True
                self.intervention_time[i] = t
                if v[i+self.ngroups]+v[i+2*self.ngroups] > maxir:
                    maxir = v[i+self.ngroups]+v[i+2*self.ngroups]
            self.maxir = maxir
        
    def subgroup_to_pop_threshold(self, v,t):
        '''
            intervention kicks in as soon as the total cumulative number of infected reaches threshold
        '''
        
        cumulativefraction = v[self.ngroups:2*self.ngroups]+v[2*self.ngroups:] 
        if np.max(cumulativefraction) > self.threshold[0] and self.int_happ[0] == False:
            for i in range(self.ngroups):
                self.int_happ[i] = True
                self.intervention_time[i] = t
    def check_time(self, v,t):
        '''
         routine that keeps track of whether and when intervention has happened
        '''
        for i in range(self.ngroups):
            if t>self.threshold[i] and t <self.threshold[i]+self.duration[i] and self.int_happ[i] == False:
                self.int_happ[i] = True
                self.intervention_time[i] = self.threshold[i]               
    
    def sir_interventionode(self,incond,t):
        '''
            routine to solve the ode
        '''
        
        sol = np.zeros((len(t)+1,len(incond)))
        sol[0,:] = incond
        dt = t[1]-t[0]
        for index,time in enumerate(t):
            
            self.check(sol[index],time)
            
            betaeff = self.beta_value(time)

        
            matrixmultiply = np.dot(betaeff, sol[index,self.ngroups: 2*self.ngroups])
            if self.ngroups==1:
                matrixmultiply= np.array([betaeff[0]*sol[index,1]])
            
            
            
            for i in range(self.ngroups):
                sol[index+1,i]  = sol[index,i]  - dt*sol[index,i]*matrixmultiply[i]
                sol[index+1,i+self.ngroups]= sol[index,i+self.ngroups]  + dt*(sol[index,i]*matrixmultiply[i] - self.gamma[i] * sol[index,i+self.ngroups]) 
                sol[index+1,i+2*self.ngroups]= sol[index,i+2*self.ngroups] +dt*(self.gamma[i]*sol[index,i+self.ngroups])
        return sol
        
    
    
    
    
    
    
    def sir(self,v,t):
        '''  
        
             DEPRECATED!!
             
             
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
        
        auxv = v.reshape(3,self.ngroups).T
        vdot = np.empty_like(auxv)
        print(auxv)
        vdot[:,0] = -np.dot(v[:,0],self.betaij)* (-v[:,1])
        vdot[:,1] = np.dot(v[:,0],self.betaij)* (-v[:,1])    - self.gamma*v[:,0]
        vdot[:2] = + self.gamma*v[:,1]
        
        vdot= vdot.T
        vdot= np.squeeze(np.asarray(vdot))
        '''
        vdot = np.empty_like(v)
        
        matrixmultiply = np.dot(self.betaij, v[self.ngroups: 2*self.ngroups])
        for i in range(self.ngroups):
            vdot[i]  = - v[i]*matrixmultiply[i]
            vdot[i+self.ngroups]= v[i]*matrixmultiply[i] - self.gamma * v[i+self.ngroups] 
            vdot[i+2*self.ngroups]= self.gamma*v[i+self.ngroups]
        return vdot  
    
    
    def sir_odeint(self, nt = 1000):
            
            t = np.linspace(0,self.tauf,nt)
            #args=(ak,dk)
            y= odeint(self.sir, self.incond, t)
            return y
    
    
    def normalepidemics(self,epsilon):
        '''
            does a normal SIR epidemic to test
            epsilon: float, initial condition in group 0
        '''
        
        Sin = np.ones(self.ngroups)
        Sin[0] = 1-epsilon
        Iin = np.zeros(self.ngroups)  
        Iin[0] = epsilon
        Rin = np.zeros(self.ngroups)
        
        self.incond = np.hstack((Sin,Iin,Rin))
        y = self.sir_odeint()
        return y




if __name__ == '__main__':
    
    
    
    #tests: with THRESHOLD

    # rate of infection matrix saved on a csv
    betaij = np.loadtxt('mixing_baseline.txt', delimiter=',')
    ngroups = len(betaij)
    #recovery rates for each group
    gamma = np.array(ngroups*[1], dtype=float)
    #final time
    tauf = 35
    SIR = SIR_model(ngroups, gamma, tauf,betaij, betain=0, betaoff=0, seed=1)
    #intervention strenght
    c = np.array(ngroups*[0.5], dtype =float)
    #initial group to be infected
    ingroup =1
    #threshold of infection 
    threshold = np.array(ngroups*[0.1], dtype=float)
    #duration of control
    duration = ngroups*[4.0]
    time =np.array(ngroups*[8])
    y=SIR.sir_intervention( c, ingroup, time, duration, nt = 3000, epsilon=0.01, intervention='time')
    
    #y[:ngroups] is the S_1(t)... S_n(t) susceptible populations evolution,
    #y[ngroups:2*ngroups] "I(t)"
    #y[2*ngroups:] "R(t)"
    
    t = np.linspace(0,tauf,3000)
  
    fig,ax = plt.subplots(3,3, figsize=(12,12))
    ax = ax.ravel()
    #plot I(t)
    for i,sub in enumerate(ax):
        #S(t)
        #sub.plot(t, y[:,i], color='b')
        #I(t)
        sub.plot(t, y[:,i+ngroups], color='r')
        #R(t)
        #sub.plot(t, y[:,i+2*ngroups], color='g')
        #intervention
        sub.vlines(SIR.intervention_time[i], 0,np.max(y[:,i+ngroups]))
        sub.set_title("group %d" %i)
    
    '''
    #test with TIME 
     # rate of infection matrix saved on a csv
    betaij = np.loadtxt('check.csv', delimiter=',')
    ngroups = len(betaij)
    #recovery rates for each group
    gamma = np.array(ngroups*[1], dtype=float)
    #final time
    tauf = 80
    SIR = SIR_model(ngroups, gamma, tauf,betaij, betain=0, betaoff=0, seed=1)
    #intervention strenght
    c = np.array(ngroups*[0.9], dtype =float)
    #initial group to be infected
    ingroup =1
    #TIME OF INFECTION
    threshold = np.array(ngroups*[4], dtype=float)
    #duration of control
    duration = ngroups*[10]
    
    y=SIR.sir_intervention( c, ingroup, threshold, duration, nt = 3000, epsilon=0.001, intervention = 'time')
    t = np.linspace(0,tauf,3000)
  
    fig,ax = plt.subplots(3,3, figsize=(12,12))
    ax = ax.ravel()
    #plot I(t)
    for i,sub in enumerate(ax):
        #S(t)
        #sub.plot(t, y[:,i], color='b')
        #I(t)
        sub.plot(t, y[:,i+ngroups], color='r')
        #R(t)
        #sub.plot(t, y[:,i+2*ngroups], color='g')
        #intervention
        sub.vlines(SIR.intervention_time[i], 0,np.max(y[:,i+ngroups]))
        sub.set_title("group %d" %i)
       
    
    '''
    
    
