import numpy as np
import numpy.random as rng
from matplotlib import rc
import matplotlib.pyplot as plt

"""
Make a plot of \pi(L1, L2)
Based on a similar plot from ABCNS
"""

rng.seed(1)
rc("font", size=18, family="serif", serif="Computer Sans")
rc("text", usetex=True)

# Resolution
N = 256
[x, y] = np.meshgrid(np.linspace(0., 5., N), np.linspace(5., 0., N))

f = np.exp(-0.5*(x-1.5)**2/1.**2)*np.exp(-0.5*((y - 5*(x/5)**2)**2)/1.**2)
f /= f.sum()

# Generate samples
M = 41
xx = 1.5 + 1.*rng.randn(M)
yy = 5*(xx/5)**2 + 1.*rng.randn(M)
keep = (xx > 0.) & (xx < 5.) & (yy > 0.) & (yy < 5.)
xx = xx[keep]
yy = yy[keep]

plt.imshow(f, extent=[x.min(), x.max(), y.min(), y.max()], cmap='Blues')
plt.plot(xx, yy, 'ko')
plt.xlabel(r'$L_1$')
plt.ylabel(r'$L_2$')
plt.title(r'Prior $\pi(L_1, L_2)$')
plt.fill_between(x[0, :][x[0, :] > 2.], 1.5, 5., color=[0, 0, 0], alpha=0.1)
plt.savefig('joint1.pdf', bbox_inches='tight')
plt.show()


ucc = np.zeros((N, N))
for i in range(N):
	for j in range(N):
		ucc[i, j] = 0
		for k in range(len(xx)):
			if (xx[k] > x[i, j]) and (yy[k] > y[i, j]):
				ucc[i, j] += 1
Xhat = ucc/len(xx)

particle_uccs = np.zeros(len(xx))
for i in range(len(xx)):
	particle_uccs[i] = 0
	for k in range(len(xx)):
		if (xx[k] >= xx[i]) and (yy[k] >= yy[i]):
			particle_uccs[i] += 1

boundary = particle_uccs == np.median(particle_uccs)
interior = particle_uccs > np.median(particle_uccs)
exterior = particle_uccs < np.median(particle_uccs)

plt.imshow(f, extent=[x.min(), x.max(), y.min(), y.max()], cmap='Blues')
plt.imshow(-(Xhat>0.25), extent=[x.min(), x.max(), y.min(), y.max()], cmap='gray', alpha=0.2)
plt.plot(xx[exterior], yy[exterior], 'go', markersize=6, label='Exterior')
plt.plot(xx[boundary], yy[boundary], 'y*', markersize=8, label='Boundary')
plt.plot(xx[interior], yy[interior], 'rs', markersize=6, label='Interior')
plt.legend(loc='upper left', numpoints=1)
plt.xlabel(r'$L_1$')
plt.ylabel(r'$L_2$')
plt.title(r'Prior $\pi(L_1, L_2)$')
plt.savefig('joint2.pdf', bbox_inches='tight')
plt.show()

