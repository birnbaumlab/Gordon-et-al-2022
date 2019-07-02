
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

def get_pprint(silent = True):
    """ Returns a silent/not silent version of print depending on original passed state """
    if silent == True:
        def pprint(*args): pass
        return pprint

    else:
        def pprint(*args): print(*args) 
        return pprint
    

