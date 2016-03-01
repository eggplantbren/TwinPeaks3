from pylab import *
import scipy

rc("font", size=18, family="serif", serif="Computer Sans")
rc("text", usetex=True)

def logsumexp(values):
	biggest = np.max(values)
	x = values - biggest
	result = np.log(np.sum(np.exp(x))) + biggest
	return result


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
  var1 = sum(WW*(scalars[:,0] - exp1)**2)
  var2 = sum(WW*(scalars[:,1] - exp2)**2)

  return [logZ, H, logWW, ess, exp1, exp2, var1, var2]

[logZ, H, logWW, ess, exp1, exp2, var1, var2] = canonical_properties(T1, T2)

#print('log(Z) = {logZ}'.format(logZ=logZ))
#print('H = {H} nats'.format(H=H))

plot(exp(logWW))
ylabel('Weight wrt canonical distribution')
title('ESS (for purposes of normalising constant calc) = {ess}'.format(ess=ess))
show()

# Resample to uniform weight
N = int(ess)
posterior_sample = zeros((N, scalars.shape[1]))
w = exp(logWW)/max(exp(logWW))
for i in range(0, N):
  while True:
    which = randint(scalars.shape[0])
    if rand() <= w[which]:
      break
  posterior_sample[i,:] = scalars[which,:]

plot(scalars[:,0], scalars[:,1], 'y.', markersize=1, label='Discarded points', alpha=1.)
hold(True)
plot(posterior_sample[:,0], posterior_sample[:,1], 'k.', markersize=1,\
					alpha=1, label='Canonical distribution')
# Adaptive axes
ss1 = sort(scalars[:,0])
ss2 = sort(scalars[:,1])
xlim([ss1[int(0.05*len(ss1))], ss1.max()])
ylim([ss2[int(0.05*len(ss2))], ss2.max()])
xlabel("$L_1$")
ylabel("$L_2$")
legend(loc='lower left', markerscale=10, numpoints=1)
savefig('output.png', bbox_inches='tight', dpi=300)
savetxt("posterior_sample.txt", posterior_sample)
show()

# Get partition function and entropy
N = 51
T1 = 10.**linspace(-2., 4., N)
T2 = T1.copy()
[T1, T2] = meshgrid(T1, T2)
T2 = T2[::-1, :]

logZ = zeros((N, N))
H = zeros((N, N))
exp1 = zeros((N, N))
exp2 = zeros((N, N))
var1 = zeros((N, N))
var2 = zeros((N, N))
depth = output[:,0].max() - output[:,0].min()

for i in range(0, N):
  for j in range(0, N):
    [logZ[i,j], H[i,j], temp1, temp2, exp1[i, j], exp2[i, j],\
				var1[i,j], var2[i,j]] = canonical_properties(T1[i, j], T2[i, j])
    # Blank out 'unreliable' results
    if H[i, j] > 0.8*depth:
      H[i, j] = NaN
      logZ[i, j] = NaN
      exp1[i, j] = NaN
      exp2[i, j] = NaN
      var1[i, j] = NaN
      var2[i, j] = NaN
  print(i+1)

figure(figsize=(11, 11))
subplot(2, 2, 1)
imshow(logZ, extent=(T1.min(), T1.max(), T2.min(), T2.max()), interpolation='nearest', cmap='viridis')
xscale('log')
yscale('log')
title(r'$\log(Z)$')
ylabel(r'$\log_{10}(T_2)$')

subplot(2, 2, 2)
imshow(H, extent=(T1.min(), T1.max(), T2.min(), T2.max()), interpolation='nearest', cmap='viridis')
xscale('log')
yscale('log')
title(r'$H$')

subplot(2, 2, 3)
imshow(exp1, extent=(T1.min(), T1.max(), T2.min(), T2.max()), interpolation='nearest', cmap='viridis')
xscale('log')
yscale('log')
title(r'$\left<\log L_1\right>$')
xlabel(r'$T_1$')
ylabel(r'$T_2$')

subplot(2, 2, 4)
imshow(exp2, extent=(T1.min(), T1.max(), T2.min(), T2.max()), interpolation='nearest', cmap='viridis')
xscale('log')
yscale('log')
title(r'$\left<\log L_2\right>$')
xlabel(r'$T_1$')

savefig('results.pdf', bbox_inches='tight')
show()

true_logZ = loadtxt('true_logZ.txt')
true_H = loadtxt('true_H.txt')
subplot(1,2,1)
err1 = logZ - true_logZ
good = logical_not(isnan(err1))
imshow(err1, extent=(log10(T1.min()), log10(T1.max()), log10(T2.min()), log10(T2.max())), interpolation='nearest', vmin=-(abs(err1[good])).max(), vmax=(abs(err1[good])).max(), cmap='coolwarm')
print("H_max = " + str(H[good].max()))
print(err1[good].min(), err1[good].max())
print(sqrt(mean(err1[good]**2)))
print(4.*H[good].max()/(abs(err1[good]).max())**2)

subplot(1,2,2)
err2 = (H - true_H)#/true_H
good = logical_not(isnan(err1))
imshow(err2, extent=(log10(T1.min()), log10(T1.max()), log10(T2.min()), log10(T2.max())), interpolation='nearest', vmin=-(abs(err2[good])).max(), vmax=(abs(err2[good])).max(), cmap='coolwarm')
print(err2[good].min(), err2[good].max())
show()

