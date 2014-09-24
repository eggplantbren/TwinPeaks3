import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
import copy

def randh(shape=None):
  if shape is not None:
    return 10.**(1.5 - 6.*rng.rand(shape))*rng.randn(shape)
  return 10.**(1.5 - 6.*rng.rand())*rng.randn()

class Model:
  N = 100
  def __init__(self):
    self.x = np.empty(Model.N)

  def from_prior(self):
    self.x = rng.rand(Model.N)
    self.calculate_scalars()

  def perturb(self):
    which = rng.randint(Model.N)
    self.x[which] += randh()
    self.x[which] = np.mod(self.x[which], 1.)
    self.calculate_scalars()
    return 0.

  def calculate_scalars(self):
    self.scalars = np.empty(2);
    self.scalars[0] = -((self.x - 0.5)**2).sum()
    self.scalars[1] = -(np.sin(4.*np.pi*self.x)**2).sum()



class Walker:
  def __init__(self):
    self.model = Model()
    self.model.from_prior()
    self.logX = rng.rand()
    self.boundary = np.array([-1E300, -1E300])

  def update(self, mcmc_steps=1000):
    self.boundary[0] = self.model.scalars[0]

    accept = 0
    for i in range(0, mcmc_steps):
      proposal = copy.deepcopy(self.model)
      logH = proposal.perturb()
      if logH > 0.:
        logH = 0.
      if rng.rand() <= np.exp(logH) and\
			(proposal.scalars[0] > self.boundary[0]):
        self.model = proposal
        accept += 1
    self.logX += np.log(rng.rand())
    print('Accepted {a}/{b}'.format(a=accept, b=mcmc_steps))

if __name__ == '__main__':
  walker = Walker()

  plt.ion()
  plt.hold(True)
  for i in range(0, 100):
    walker.update()
    plt.plot(walker.logX, walker.boundary[0], 'bo')
    plt.draw()

  plt.ioff()
  plt.show()

