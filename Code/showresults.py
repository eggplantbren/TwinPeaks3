from pylab import *
from postprocess import logsumexp

scalars = loadtxt('scalars.txt')
logw = loadtxt('logw.txt')
smallest = min([scalars.shape[0], logw.size])
scalars = scalars[0:smallest, :]
logw = logw[0:smallest]

# Prior weights, normalised
logw = logw - logsumexp(logw)

# Posterior weights, unnormalised
logW = logw + 10*scalars[:,0] + scalars[:,1]

# Normaliser
logZ = logsumexp(logW)
print('log(Z) = {logZ}'.format(logZ=logZ))

# Posterior weights, normalised
logWW = logW - logZ
ess = -sum(exp(logWW)*logWW)

plot(exp(logWW))
ylabel('Posterior Weight')
title('ESS = {ess}'.format(ess=ess))
show()

