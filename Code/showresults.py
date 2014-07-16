from pylab import *
from postprocess import logsumexp

# First calculate things about the scalars (e.g. the normalising constant)
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

# Posterior weights, normalised
logWW = logW - logZ
ess = exp(-sum(exp(logWW)*logWW))

# Information
H = sum(exp(logWW)*(logWW - logw))

print('log(Z) = {logZ}'.format(logZ=logZ))
print('H = {H} nats'.format(H=H))

plot(exp(logWW))
ylabel('Weight wrt canonical distribution')
title('ESS (for purposes of normalising constant calc) = {ess}'.format(ess=ess))
show()

## Now do it all again with the thinned samples
#scalars = loadtxt('scalars_thinned.txt')
#logw = loadtxt('logw_thinned.txt')
#smallest = min([scalars.shape[0], logw.size])
#scalars = scalars[0:smallest, :]
#logw = logw[0:smallest]

## Prior weights, normalised
#logw = logw - logsumexp(logw)

## Posterior weights, unnormalised
#logW = logw + 10*scalars[:,0] + scalars[:,1]

## Normaliser
#logZ = logsumexp(logW)

## Posterior weights, normalised
#logWW = logW - logZ
#ess = exp(-sum(exp(logWW)*logWW))

