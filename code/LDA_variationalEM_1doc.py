# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 19:47:01 2015

@author: marcus

Test of variational EM algorithm for LDA with 1 document.
"""

import numpy as np
import numpy.random as rng
from scipy.special import psi as digamma

W = 4 #number of "words" (= bins)
K = 2 #number of topics


B0 = np.ones((K,W), dtype=float)
B0[0] = 1000 * B0[0]  # we have lots of confidence about topic 0: it's the background...

alpha = np.ones((K,1), dtype=float) # hyperparams for Dir giving mixture coeffs
gamma = np.ones((K,1), dtype=float) # initial params for the Dirichlet topic distribution
beta = rng.random(size=(K,W)) # initial parameters of Multinomial word distribution, for each topic
beta = beta / beta.sum(1).reshape(K,1)   # normalised along the rows
phi = np.ones((K,W), dtype=float)

# SOME MADE-UP DATA........................
# EG: 50 pixels from background, and another 150 from [.5, .5, 0, 0]
# So mixing proportion is about 1/4 background.
n = np.array([205, 195, 50, 50]).reshape((1,W))  # the data as "word" counts

# Run variational EM algorithm to find the best-guess gammas.
for t in range(100):
    # E step: update phi (normalised) and gamma
    for tt in range(10):
        phi = beta * np.exp(digamma(gamma) - digamma(gamma.sum()))
        phi = phi / phi.sum(1).reshape(K,1)   # normalised along the rows
        gamma = alpha + (n * phi).sum(1).reshape(K,1)
    
    # M step: update beta (normalised)
    beta = B0 + (n * phi)
    beta =  1. / beta.sum(1).reshape(K,1) * beta   # normalised along the rows

print 'phi: '
print phi
print 'gamma:'
print gamma
print 'most likely proportion of Background: ', gamma[0]/gamma.sum()
print 'beta :'
print beta

# Observations: both gamma and beta are wrong!
# Mix proportion 0.38, should be 0.25
# The beta distribution for Source is really skewed by the slight difference in 205 versus 195.
# That beta skew is really weird: should step it through on paper maybe.
#
# Also, Blei et al different from Murphy in the second digamma term?? But I've tried both and it 
# doesn't solve the problem anyway. :-(