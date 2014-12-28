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

  return [logZ, H, logWW, ess, exp1, exp2]

[logZ, H, logWW, ess, exp1, exp2] = canonical_properties(T1, T2)

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

# Get partition function and entropy
N = 51
T1 = 10.**linspace(-2., 5., N)
T2 = T1.copy()
[T1, T2] = meshgrid(T1, T2)
T2 = T2[::-1, :]

logZ = zeros((N, N))
H = zeros((N, N))
exp1 = zeros((N, N))
exp2 = zeros((N, N))
depth = output[:,0].max() - output[:,0].min()

for i in xrange(0, N):
  for j in xrange(0, N):
    [logZ[i,j], H[i,j], temp1, temp2, exp1[i, j], exp2[i, j]] = canonical_properties(T1[i, j], T2[i, j])
    # Blank out 'unreliable' results
    if H[i, j] > depth/1.5:
      H[i, j] = NaN
      logZ[i, j] = NaN
      exp1[i, j] = NaN
      exp2[i, j] = NaN
  print(i+1)

subplot(2, 2, 1)
imshow(logZ, extent=(T1.min(), T1.max(), T2.min(), T2.max()), interpolation='nearest')
title('log(Z)')

subplot(2, 2, 2)
imshow(H, extent=(T1.min(), T1.max(), T2.min(), T2.max()), interpolation='nearest')
title('H')

subplot(2, 2, 3)
imshow(exp1, extent=(T1.min(), T1.max(), T2.min(), T2.max()), interpolation='nearest')
title('<S1>')

subplot(2, 2, 4)
imshow(exp2, extent=(T1.min(), T1.max(), T2.min(), T2.max()), interpolation='nearest')
title('<S2>')

show()

