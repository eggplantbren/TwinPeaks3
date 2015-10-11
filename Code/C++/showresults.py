from pylab import *
import scipy
from postprocess import logsumexp

# Temperatures
T1, T2 = 0.1, 1.

# First calculate things about the scalars (e.g. the normalising constant)
output = loadtxt('sample_info.txt')
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

plot(posterior_sample[:,0], posterior_sample[:,1], 'b.', markersize=5, alpha=0.5)
hold(True)
plot(scalars[:,0], scalars[:,1], 'r.', markersize=1, alpha=0.2)
xlabel("Scalar 1")
ylabel("Scalar 2")
savetxt("posterior_sample.txt", posterior_sample)
show()

# Get partition function and entropy
N = 51
T1 = 10.**linspace(-2.5, 5., N)
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
    if H[i, j] > 0.7*depth:
      H[i, j] = NaN
      logZ[i, j] = NaN
      exp1[i, j] = NaN
      exp2[i, j] = NaN
  print(i+1)

subplot(2, 2, 1)
imshow(logZ, extent=(log10(T1.min()), log10(T1.max()), log10(T2.min()), log10(T2.max())), interpolation='nearest', cmap='Purples')
title('log(Z)')

subplot(2, 2, 2)
imshow(H, extent=(log10(T1.min()), log10(T1.max()), log10(T2.min()), log10(T2.max())), interpolation='nearest', cmap='Purples')
title('H')

subplot(2, 2, 3)
imshow(exp1, extent=(log10(T1.min()), log10(T1.max()), log10(T2.min()), log10(T2.max())), interpolation='nearest', cmap='Purples')
title('<S1>')

subplot(2, 2, 4)
imshow(exp2, extent=(log10(T1.min()), log10(T1.max()), log10(T2.min()), log10(T2.max())), interpolation='nearest', cmap='Purples')
title('<S2>')

show()

true_logZ = loadtxt('true_logZ.txt')
true_H = loadtxt('true_H.txt')
subplot(1,2,1)
err1 = logZ - true_logZ
good = logical_not(isnan(err1))
imshow(err1, extent=(log10(T1.min()), log10(T1.max()), log10(T2.min()), log10(T2.max())), interpolation='nearest', cmap='coolwarm')
print("H_max = " + str(H[good].max()))
print(err1[good].min(), err1[good].max())
print(4.*H[good].max()/(abs(err1[good]).max())**2)

subplot(1,2,2)
err2 = (H - true_H)#/true_H
good = logical_not(isnan(err1))
imshow(err2, extent=(log10(T1.min()), log10(T1.max()), log10(T2.min()), log10(T2.max())), interpolation='nearest', cmap='coolwarm')
print(err2[good].min(), err2[good].max())
show()


