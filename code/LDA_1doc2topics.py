# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 13:18:34 2015
@author: marcus

Testing my idea of variational EM on LDA with ONE document, 2 topics!
"""

import numpy as np
import numpy.random as rng
from scipy.special import psi as digamma

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
print 'pi: ', pi
print B

for t in range(5):
    di = digamma(pi[0])-digamma(pi[1])
    di_array = -1. * np.array([di, -di]).reshape((2,1))
    C = n * B * np.exp(di_array)
    print 'C :'
    print C

    # E step
    pi = alpha + C.sum(1).reshape((2,1))
    pi = pi / np.sum(pi)  # normalise it
    print 'pi: ', pi
    
    # M step
    B = gamma + C.sum(0).reshape((1,K))
    B =  1. / B.sum(1).reshape(2,1) * B   # ugly, but now normalised along the rows
    print 'B :'
    print B

    
