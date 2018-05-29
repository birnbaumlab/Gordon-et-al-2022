
# 

# nonstandard libraries
import numpy as np
import distance

# FACTORY METHODS

def hamming_dist_np(a,b):
    """ Hamming distance of np arrays (FAST) """
    return np.count_nonzero(a!=b)

def hamming_dist(a,b):
    """ Hamming distance """
    return distance.hamming(a,b)

def levenshtein_dist(a,b):
    """ Levenshtein distance """
    return distance.levenshtein(a,b)

