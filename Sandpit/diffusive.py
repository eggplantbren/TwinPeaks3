import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
import copy

def randh(shape=None):
	if shape is not None:
		return 10.**(1.5 - 6.*rng.rand(shape))*rng.randn(shape)
	return 10.**(1.5 - 6.*rng.rand())*rng.randn()

class Walker:
	N = 100
	def __init__(self):
		self.x = np.empty(Walker.N)

	def from_prior(self):
		self.x = rng.rand(Walker.N)
		self.calculate_scalars()

	def perturb(self):
		which = rng.randint(Walker.N)
		self.x[which] += randh()
		self.x[which] = np.mod(self.x[which], 1.)
		self.calculate_scalars()
		return 0.

	def calculate_scalars(self):
		self.scalars = np.empty(2);
		self.scalars[0] = -((self.x - 0.5)**2).sum()
		self.scalars[1] = -(np.sin(4.*np.pi*self.x)**2).sum()

	def update(self, mcmc_steps=100):
		accept = 0
		walker = copy.deepcopy(self)
		for i in range(0, mcmc_steps):
			proposal = copy.deepcopy(walker)
			logH = proposal.perturb()
			if logH > 0.:
				logH = 0.
			if rng.rand() <= np.exp(logH):
				walker = proposal
				accept += 1
		return walker

# Calculate a corner mass
def corner_mass(s, keep):
	in_corner = np.logical_and(keep[:,0] < s[0], keep[:,1] < s[1])
	return float(in_corner.sum())/(len(keep) - 1)

if __name__ == '__main__':
	walker = Walker()
	walker.from_prior()

	# Store scalars from the prior
	keep = np.zeros((100, 2))

	plt.ion()
	plt.hold(True)
	for i in range(0, 100):
		walker = walker.update()
		keep[i, :] = walker.scalars
		plt.plot(walker.scalars[0], walker.scalars[1], 'bo', alpha=0.2)
		plt.title(i+1)
		plt.draw()

	plt.plot(-8., -50., 'ro')
	# Calculate corner mass
	print(corner_mass(np.array([-8., -50.]), keep))

	plt.ioff()
	plt.show()


