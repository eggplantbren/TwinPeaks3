from pylab import *

x = linspace(0., 1., 2001)

def truth(T1, T2, do_plot=False):
  p = exp(-(x - 0.5)**2/T1 - sin(2.*pi*x/0.5)**2/T2)
  Z = trapz(p, x=x)
  H = trapz(p/Z*log(p/Z), x=x)

  logZ = 200*log(Z)
  H *= 200

  if do_plot:
    plot(x, p)
    title([logZ, H])
    show()

  return [logZ, H]
  

def grid():
  # Calculate log(Z) and H for some canonical distributions
  T1 = 10.**(linspace(-2., 5., 51))
  T2 = T1.copy()
  [T1, T2] = meshgrid(T1, T2)
  T2 = T2[::-1, :]
  logZ = T1.copy()
  H = T1.copy()

  for i in xrange(0, 51):
    for j in xrange(0, 51):
      [logZ[i, j], H[i, j]] = truth(T1[i, j], T2[i, j])
    print(i+1)

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
  title('H, max={m}'.format(m=H.max()))

  show()
  return [logZ, H]

if __name__ == '__main__':
  truth(0.1, 1., do_plot=True)
  [logZ, H] = grid()

  savetxt('true_logZ.txt', logZ)
  savetxt('true_H.txt', H)


