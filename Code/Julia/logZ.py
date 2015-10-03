"""
Load output.txt and compute logZ
"""

import numpy as np
import matplotlib.pyplot as plt
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

temperatures = [0.1, 1.]

print('log(Z) = {a}'.format(a=logZ(temperatures)))
print('H = {h}'.format(h=H(temperatures)))

logp = output[:,0] + output[:,1]/temperatures[0] + output[:,2]/temperatures[1]\
				- logZ(temperatures)
plt.plot(np.exp(logp - logp.max()))
plt.show()
