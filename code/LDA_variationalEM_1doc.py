# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 19:47:01 2015

@author: marcus

Test of variational EM algorithm for LDA with 1 document.
"""

import numpy as np
#import numpy.random as rng
from scipy.special import psi as digammaFn
np.set_printoptions(precision=3)

V = 4 #number of "words" in the vocabulary (= bins)
K = 2 #number of topics


beta = np.ones((K,V), dtype=float)
beta[0] = 1000 * beta[0]  # we have lots of confidence about topic 0: it's the background...

alpha = np.ones((K,1), dtype=float) # hyperparams for Dir giving mixture coeffs
alphaTilde = np.ones((K,1), dtype=float) # initial params for the Dirichlet topic distribution
B = np.ones((K,V), dtype=float) # initial parameters of Multinomial word distribution, for each topic
B = B / B.sum(1).reshape(K,1)   # normalised along the rows
cTilde = np.ones((K,V), dtype=float)

# SOME MADE-UP DATA........................
# EG: 50 pixels from background, and another 150 from [.5, .5, 0, 0]
# So mixing proportion is about 1/4 background.
n = np.array([205, 200, 50, 50]).reshape((1,V))  # the data as "word" counts

# Run variational EM algorithm to find the best-guess alphaTildes.
for t in range(5):
    # E step: update cTilde (normalised) and alphaTilde
    for tt in range(1):
        theFactor = np.exp(digammaFn(alphaTilde) - digammaFn(alphaTilde.sum()))
        cTilde = B * theFactor
        cTilde = cTilde / cTilde.sum(1).reshape(K,1)   # normalised along the rows
        alphaTilde = alpha + (n * cTilde).sum(1).reshape(K,1)
        print 'theFactor.shape is ',theFactor.shape,' but therefore is having NO effect!!'
        print 'so..... I have some dimensionality problem to understand here'
    
    # M step: update beta (normalised)
    B = beta + (n * cTilde)
    B =  B / B.sum(1).reshape(K,1)   # normalised along the rows

print 'cTilde: '
print cTilde
print 'alphaTilde:'
print alphaTilde
print 'most likely proportion of Background: ', alphaTilde[0]/alphaTilde.sum()
print 'B :'
print B

# Observations: both alphaTilde and B are wrong!
# Mix proportion 0.38, should be 0.25
# The B distribution for Source is really skewed by the slight difference in 205 versus 195.
# That B skew is really weird: should step it through on paper maybe.
#
# Also, Blei et al different from Murphy in the second dialphaTilde term?? But I've tried both and it 
# doesn't solve the problem anyway. :-(