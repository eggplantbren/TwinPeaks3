"""
Code to load results and resample to generate a sample
from a canonical distribution.
"""

from pylab import *
import scipy
import colormaps as cmaps

rc("font", size=18, family="serif", serif="Computer Sans")
rc("text", usetex=True)

def logsumexp(values):
	biggest = np.max(values)
	x = values - biggest
	result = np.log(np.sum(np.exp(x))) + biggest
	return result


# Temperatures
T1, T2 = 2., 10000.

# First calculate things about the scalars (e.g. the normalising constant)
sample = loadtxt('sample.txt')
scalars = sample[:,1:3]
logw = sample[:,0]
smallest = min([scalars.shape[0], logw.size])
scalars = scalars[0:smallest, :]
logw = logw[0:smallest]

# Prior weights, normalised
logw = logw - logsumexp(logw)

# Evaluate normalising constant at any temperature
def canonical_properties(temperature1, temperature2):
  # Posterior weights, unnormalised
  logW = logw + scalars[:,0]/temperature1 + scalars[:,1]/temperature2

  # Normaliser
  logZ = logsumexp(logW)

  # Posterior weights, normalised
  logWW = logW - logZ
  WW = exp(logWW)
  ess = exp(-sum(WW*logWW))

  # Information
  H = sum(WW*(logWW - logw))

  # Expected values of scalars
  exp1 = sum(WW*scalars[:,0])
  exp2 = sum(WW*scalars[:,1])
  var1 = sum(WW*(scalars[:,0] - exp1)**2)
  var2 = sum(WW*(scalars[:,1] - exp2)**2)

  return [logZ, H, logWW, ess, exp1, exp2, var1, var2]

[logZ, H, logWW, ess, exp1, exp2, var1, var2] = canonical_properties(T1, T2)
WW = exp(logWW - logWW.max())

plot(WW)
ylabel('Weight wrt canonical distribution')
title('ESS (for purposes of normalising constant calc) = {ess}'.format(ess=ess))
show()

# Resample to uniform weight
N = int(ess)
posterior_sample = zeros((N, sample.shape[1]))

for i in range(0, N):
	while True:
		which = randint(sample.shape[0])
		if rand() <= WW[which]:
			break
	posterior_sample[i,:] = sample[which,:]

savetxt("canonical_sample.txt", posterior_sample)

