import numpy as np
import numpy.random as rng
from matplotlib import rc
import matplotlib.pyplot as plt

"""
Make a plot of \pi(L1, L2)
Based on a similar plot from ABCNS
"""

rng.seed(0)
rc("font", size=18, family="serif", serif="Computer Sans")
rc("text", usetex=True)

# Resolution
N = 256
[x, y] = np.meshgrid(np.linspace(0., 5., N), np.linspace(5., 0., N))

f = np.exp(-0.5*(x-3.5)**2/1.**2)*np.exp(-0.5*((y - 5*(x/5)**2)**2)/0.3**2)
f /= f.sum()

# Generate samples
M = 20
xx = 3.5 + 1.*rng.randn(M)
yy = 5*(xx/5)**2 + 0.3*rng.randn(M)
keep = (xx > 0.) & (xx < 5.) & (yy > 0.) & (yy < 5.)
xx = xx[keep]
yy = yy[keep]

plt.imshow(f, extent=[x.min(), x.max(), y.min(), y.max()], cmap='Blues')
plt.plot(xx, yy, 'ko')
plt.xlabel(r'$L_1$')
plt.ylabel(r'$L_2$')
plt.title(r'Prior $\pi(L_1, L_2)$')
plt.axhline(1.3095, xmin=0., xmax=2.6455/5., color='k')
plt.axvline(2.6455, ymin=0., ymax=1.3095/5., color='k')
plt.axhline(1.3095, linestyle='--', color='k')
plt.axvline(2.6455, linestyle='--', color='k')
plt.axhline(0.841, xmin=0., xmax=2.6455/5., linestyle='-.', color='k')
plt.fill_between(x[0, :][x[0, :] < 2.6455], 0., 1.3095, color=[0.6, 0.6, 0.6], alpha=0.2)
plt.savefig('joint.pdf', bbox_inches='tight')
plt.show()

