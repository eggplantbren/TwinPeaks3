from pylab import *
from postprocess import logsumexp

scalars = loadtxt('scalars.txt')
logw = loadtxt('logw.txt')

# Prior weights, normalised
logw = logw - logsumexp(logw)

# Posterior weights, unnormalised
logW = logw + scalars[:,0]

# Normaliser
logZ = logsumexp(logW)
print('log(Z) = {logZ}'.format(logZ=logZ))

# Posterior weights, normalised
logWW = logW - logZ

plot(exp(logWW))
ylabel('Posterior Weight')
show()

