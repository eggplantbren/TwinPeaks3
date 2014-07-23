from pylab import *

x = linspace(0., 1., 2001)

def truth(L1, L2, do_plot=False):
  p = exp(-L1*(x - 0.5)**2 - L2*sin(2.*pi*x/0.5)**2)
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
  # The two lagrange multipliers
  L1 = linspace(0., 20., 101)
  L2 = linspace(0., 20., 101)
  [L1, L2] = meshgrid(L1, L2)
  L2 = L2[::-1, :]

  logZ = L1.copy()
  info = L1.copy()

  for i in xrange(0, 101):
    for j in xrange(0, 101):
      [logZ[i, j], info[i, j]] = truth(L1[i, j], L2[i, j])
    print(i+1)

  subplot(1, 2, 1)
  imshow(logZ)

  subplot(1, 2, 2)
  imshow(info)
  show()


if __name__ == '__main__':
  truth(10., 1., do_plot=True)
  grid()


