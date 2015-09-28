"""
Various useful functions
"""

import numpy as np

def logsumexp(values):
	biggest = np.max(values)
	x = values - biggest
	result = np.log(np.sum(np.exp(x))) + biggest
	return result

def logdiffexp(x1, x2):
	biggest = x1
	xx1 = x1 - biggest
	xx2 = x2 - biggest
	result = np.log(np.exp(xx1) - np.exp(xx2)) + biggest
	return result

