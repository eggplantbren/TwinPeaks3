"""
Load sample_info.txt and compute logZ
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

temperatures = [0.1, 1.]
print('log(Z) = {a}'.format(a=logZ(temperatures)))
print('H = {h}'.format(h=H(temperatures)))

import matplotlib.pyplot as plt
plt.figure(1)
plt.plot(output[:,1], output[:,2], 'b.', markersize=1)
plt.xlabel('Scalar 1')
plt.ylabel('Scalar 2')

plt.figure(2)
logweight = output[:,0] + output[:,1]/temperatures[0] + output[:,2]/temperatures[1]
plt.plot(np.exp(logweight - logweight.max()))
plt.xlabel('Iteration')
plt.ylabel('Relative posterior weight')
plt.show()

