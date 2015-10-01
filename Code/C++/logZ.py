"""
Load output.txt and compute logZ
"""

import numpy as np
from Utils import logsumexp

# Load the data
output = np.loadtxt('sample_info.txt')

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

print('log(Z) = {a}'.format(a=logZ([0.1, 1.])))
print('H = {h}'.format(h=H([0.1, 1.])))

