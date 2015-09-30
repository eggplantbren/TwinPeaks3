"""
Load output.txt and compute logZ
"""

import numpy as np
from Utils import logsumexp

# Load the data
output = np.loadtxt('output.txt')

def logZ(temperatures):
	"""
	Compute the normalising constant. Only works for two scalars.
	"""
	return logsumexp(output[:,0] + output[:,1]/temperatures[0] + output[:,2]/temperatures[1])

def H(temperatures):
	"""
	Compute the information. Only works for two scalars.
	"""
	# Log posterior weights (normalised)
	logp = output[:,0] + output[:,1]/temperatures[0] + output[:,2]/temperatures[1]\
				- logZ(temperatures)
	return np.sum(np.exp(logp)*(logp - output[:,0]))

print(logZ([0.1, 1.]))
print(H([0.1, 1.]))

