import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
from copy import deepcopy

from Walker import Walker

class Sampler:
	"""
	An object of this class is a sampler.
	"""

	def __init__(self, num_particles):
		"""
		Constructor: pass in the number of particles
		"""
		self.num_particles = num_particles
		self.walkers = [Walker() for i in range(0, num_particles)]

	def initialise(self):
		"""
		Generate all the walkers from the prior
		"""
		self.all_scalars = []	# Values of scalars for all walkers
		for walker in self.walkers:
			walker.from_prior()
			self.all_scalars.append(walker.scalars)
		self.all_scalars = np.array(self.all_scalars)
		return


sampler = Sampler(1000)
sampler.initialise()
print(sampler.all_scalars.shape)

