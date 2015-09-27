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

		# Open and close output file to clear it
		f = open('output.txt', 'w')
		f.close()

	def initialise(self):
		"""
		Generate all the walkers from the prior
		"""
		self.all_scalars = []	# Values of scalars for all walkers
		for walker in self.walkers:
			walker.from_prior()
			self.all_scalars.append(walker.scalars)
		self.all_scalars = np.array(self.all_scalars)
		self.iteration = 0
		self.log_prior_mass = 0. # Fraction of prior mass still in play
		return

	def do_iteration(self):
		"""
		Do an NS iteration
		"""
		counts = self.rectangle_counts

		# Choose one with minimum count to discard
		temp = np.nonzero(counts == counts.min())[0]
		which = temp[rng.randint(len(temp))]

		# Estimate the fraction of remaining prior mass being eliminated
		frac = (1. + counts[which])/self.num_particles

		# Write discarded particle info to disk
		# log prior mass estimate, then scalars
		f = open('output.txt', 'a')
		logw = np.log(float(self.iteration)/self.num_particles) + \
				self.iteration*np.log(1. - np.log(float(self.iteration)/self.num_particles))
		line = str(logw)
		for s in self.all_scalars[which, :]:
			line += str(s) + ' '
		f.write(line + '\n')
		f.close()

	@property
	def rectangle_counts(self):
		"""
		Count how many walkers are within the rectangle
		of each walker.
		"""
		counts = np.empty(self.num_particles, dtype='int64')
		for i in range(0, self.num_particles):
			counts[i] = np.sum(Sampler.is_in_rectangle\
									(self.all_scalars, self.all_scalars[i]))
		return counts

	@staticmethod
	def is_in_rectangle(scalars, rectangle):
		"""
		Rectangle should be a numpy array of length 2
		Scalars can be of shape (a, 2)
		Returns true for each row of 'scalars' that is within the rectangle
		defined by 'rectangle'.
		"""
		return np.all(scalars < rectangle, axis=1)

sampler = Sampler(1000)
sampler.initialise()
sampler.do_iteration()


