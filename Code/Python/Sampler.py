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
		for walker in Walkers:
			walker.from_prior()
		return


