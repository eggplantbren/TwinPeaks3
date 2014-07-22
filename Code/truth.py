from pylab import *

def truth():
  x = linspace(0., 1., 2001)

  # The two lagrange multipliers
  L1 = linspace(0., 20., 101)
  L2 = linspace(0., 20., 101)
  [L1, L2] = meshgrid(L1, L2)
  L2 = L2[::-1, :]

  logZ = L1.copy()
  info = L1.copy()

  for i in xrange(0, 101):
    for j in xrange(0, 101):
      p = exp(-L1[i, j]*(x - 0.5)**2 - L2[i, j]*sin(2.*pi*x/0.5)**2)
      Z = trapz(p, x=x)
      H = trapz(p/Z*log(p/Z), x=x)

      logZ[i, j] = 200*log(Z)
      info[i, j] = H
    print(i+1)

  subplot(1, 2, 1)
  imshow(logZ)

  subplot(1, 2, 2)
  imshow(info)
  show()


if __name__ == '__main__':
  truth()

