from pylab import *
import colormaps as cmaps

rc("font", size=18, family="serif", serif="Computer Sans")
rc("text", usetex=True)

x = linspace(0., 1., 10001)

def truth(T1, T2, do_plot=False):
  p = exp(-(x - 0.5)**2/T1 - x/T2)#sin(4.*pi*x)**2/T2)
  Z = trapz(p, x=x)
  H = trapz(p/Z*log(p/Z + 1E-300), x=x)

  logZ = 100*log(Z)
  H *= 100

  if do_plot:
    plot(x, p)
    title([logZ, H])
    show()

  return [logZ, H]
  

def grid():
  # Calculate log(Z) and H for some canonical distributions
  T1 = 10.**(linspace(-2., 4., 51))
  T2 = T1.copy()
  [T1, T2] = meshgrid(T1, T2)
  T2 = T2[::-1, :]
  logZ = T1.copy()
  H = T1.copy()

  for i in range(0, 51):
    for j in range(0, 51):
      [logZ[i, j], H[i, j]] = truth(T1[i, j], T2[i, j])
    print(i+1)

  figure(1)
  subplot(1, 2, 1)
  imshow(logZ, extent=(T1.min(), T1.max(), T2.min(), T2.max()), interpolation='nearest', cmap=cmaps.viridis)
  xlabel(r'$T_1$')
  ylabel(r'$T_2$')
  title(r'$\ln(Z)$')
  xscale('log')
  yscale('log')

  subplot(1, 2, 2)
  imshow(H, extent=(T1.min(), T1.max(), T2.min(), T2.max()), interpolation='nearest', cmap=cmaps.viridis)
  xlabel(r'$T_1$')
  title(r'$H$')
  xscale('log')
  yscale('log')
  print('H_max = ', H.max())
  savefig('truth.pdf', bbox_inches='tight')
  show()
  return [logZ, H]

if __name__ == '__main__':
  truth(0.1, 1., do_plot=True)
  [logZ, H] = grid()

  savetxt('true_logZ.txt', logZ)
  savetxt('true_H.txt', H)


