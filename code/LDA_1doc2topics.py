# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 13:18:34 2015
@author: marcus

Testing my idea of variational EM on LDA with ONE document, 2 topics!
"""

import numpy as np
import numpy.random as rng

K = 4 #number of bins
n = np.array([10, 10, 50, 50])  # the counts
gamma = np.ones((2,K), dtype=float) # hyperparams for the bin counts
gamma[0] = gamma[0]*100
alpha = np.ones((2,1), dtype=float) # hyperparams for Dir giving mixture coeffs

# initialisation of pi and B
pi = np.ones((2,1), dtype=float)
pi = pi / np.sum(pi)  # normalise it

B = np.ones((2,K), dtype=float)
B =  1. / B.sum(1).reshape(2,1) * B   # ugly, but now normalised along the rows
#print pi
#print B

for t in range(5):
    print pi
    print B
    print '\n'
    # E step
    pi = alpha + (np.transpose(pi) * np.dot(B, n)).reshape(alpha.shape)
    pi = pi / np.sum(pi)  # normalise it
    #print pi
    
    # M step
    B = gamma + n.reshape(1, K) * B * pi
    B =  1. / B.sum(1).reshape(2,1) * B   # ugly, but now normalised along the rows

    