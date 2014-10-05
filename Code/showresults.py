from pylab import *
import scipy
from postprocess import logsumexp

# Temperatures
T1, T2 = 0.3, 0.3

# First calculate things about the scalars (e.g. the normalising constant)
scalars = loadtxt('scalars.txt')
logw = loadtxt('logw.txt')
smallest = min([scalars.shape[0], logw.size])
scalars = scalars[0:smallest, :]
logw = logw[0:smallest]

# Prior weights, normalised
logw = logw - logsumexp(logw)

# Posterior weights, unnormalised
logW = logw + scalars[:,0]/T1 + scalars[:,1]/T2

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



# Calculate log(Z) and H for some canonical distributions
T1 = exp(linspace(-3., 10., 101))
T2 = exp(linspace(-3., 10., 101))
[T1, T2] = meshgrid(T1, T2)
T2 = T2[::-1, :]
logZ = T1.copy()
H = T1.copy()
ess = T1.copy()

for i in xrange(0, 101):
  for j in xrange(0, 101):
    # Prior weights, normalised
    logw = logw - logsumexp(logw)

    # Posterior weights, unnormalised
    logW = logw + scalars[:,0]/T1[i, j] + scalars[:,1]/T2[i, j]

    # Normaliser
    logZ[i, j] = logsumexp(logW)

    # Posterior weights, normalised
    logWW = logW - logZ[i, j]
    ess[i,j] = exp(-sum(exp(logWW)*logWW))

    # Information
    H[i, j] = sum(exp(logWW)*(logWW - logw))
  print(i+1)

logZ[ess < 5] = nan
H[ess < 5] = nan

figure(1)
subplot(1, 2, 1)
imshow(logZ, extent=(log(T1.min()), log(T1.max()), log(T2.min()), log(T2.max())))
xlabel('ln($T_1$)')
ylabel('ln($T_2$)')
title('log(Z)')

subplot(1, 2, 2)
imshow(H, extent=(log(T1.min()), log(T1.max()), log(T2.min()), log(T2.max())))
xlabel('ln($T_1$)')
ylabel('ln($T_2$)')
title('H, max={m}'.format(m=H[logical_not(isnan(H))].max()))

figure(2)
imshow(ess, extent=(log(T1.min()), log(T1.max()), log(T2.min()), log(T2.max())))
xlabel('ln($T_1$)')
ylabel('ln($T_2$)')
title('ESS, max = {m}'.format(m=ess.max()))

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

