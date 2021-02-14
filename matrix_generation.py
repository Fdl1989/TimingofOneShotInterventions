#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 09:41:50 2020

@author: francesco
"""

import numpy as np
import sys

if __name__ == '__main__':
    '''
     This code accepts a number in input to seed the rng (or time), 
     then it generates a random matrix with rules specified in the comments
    '''
    if len(sys.argv)>1:
        seed = int(sys.argv[1])
        np.random.seed(seed)
    else:
        import time 
        seed = time.time()
        np.random.seed()
    '''
    9 populations
    '''
    no_population = 9 
    matrix = np.zeros((no_population,no_population))
    
    '''
    diagonal = 2 + unif(0,1) -0.5
    '''
    Diag_mat = 1.5 + np.random.uniform(0,1, no_population) 
    M = np.diag(Diag_mat)
    
    '''
    off_diag = max_diagonal/10 * unif(0,1) 
    '''
    max_elem = np.max(M)
    upper_diag = 0.1*max_elem*(np.random.uniform(0,1,no_population-1))
    lower_diag = 0.1*max_elem*(np.random.uniform(0,1,no_population-1))
    
    M += np.diag(upper_diag, k=1)
    M += np.diag(lower_diag, k=-1)
    
    np.savetxt("mat_realizations/%d.txt"%seed, M)