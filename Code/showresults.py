from pylab import *
import scipy
from postprocess import logsumexp

# Temperatures
T1, T2 = 0.1, 1.

# First calculate things about the scalars (e.g. the normalising constant)
output = loadtxt('output.txt')
scalars = output[:,1:]
logw = output[:,0]
smallest = min([scalars.shape[0], logw.size])
scalars = scalars[0:smallest, :]
logw = logw[0:smallest]

# Prior weights, normalised
logw = logw - logsumexp(logw)

# Evaluate normalising constant at any temperature
def logZ(temperature1, temperature2):
  # Posterior weights, unnormalised
  logW = logw + scalars[:,0]/temperature1 + scalars[:,1]/temperature2

  # Normaliser
  logZ = logsumexp(logW)

  # Posterior weights, normalised
  logWW = logW - logZ
  ess = exp(-sum(exp(logWW)*logWW))

  # Information
  H = sum(exp(logWW)*(logWW - logw))
  return [logZ, H, logWW, ess]

[logZ, H, logWW, ess] = logZ(T1, T2)

print('log(Z) = {logZ}'.format(logZ=logZ))
print('H = {H} nats'.format(H=H))

plot(exp(logWW))
ylabel('Weight wrt canonical distribution')
title('ESS (for purposes of normalising constant calc) = {ess}'.format(ess=ess))
show()

# Resample to uniform weight
N = int(ess)
posterior_sample = zeros((N, scalars.shape[1]))
w = exp(logWW)/max(exp(logWW))
for i in xrange(0, N):
  while True:
    which = randint(scalars.shape[0])
    if rand() <= w[which]:
      break
  posterior_sample[i,:] = scalars[which,:]

plot(posterior_sample[:,0], posterior_sample[:,1], 'b.', markersize=1)
savetxt("posterior_sample.txt", posterior_sample)
show()


