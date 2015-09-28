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

		# Rectangles that have been forbidden so far
		self.forbidden_rectangles = None

		# Number of iterations done
		self.iteration = 0

		# Fraction of prior mass still in play
		self.log_prior_mass = 0.

	def initialise(self):
		"""
		Generate all the walkers from the prior
		"""
		self.all_scalars = []	# Values of scalars for all walkers
		for walker in self.walkers:
			walker.from_prior()
			self.all_scalars.append(walker.scalars)
		self.all_scalars = np.array(self.all_scalars)
		self.forbidden_rectangles = np.empty((0, self.all_scalars.shape[1]))
		return

	def do_iteration(self):
		"""
		Do an NS iteration
		"""
		counts = self.rectangle_counts

		# Choose one with minimum count to discard
		temp = np.nonzero(counts == counts.min())[0]
		which = temp[rng.randint(len(temp))]

		# Need to figure out what happens in this case
		if counts.min() != 0:
			print("counts.min() != 0")

		# Estimate the fraction of *remaining* prior mass being eliminated
		frac = float(1 + counts[which])/self.num_particles

		# Write discarded particle info to disk
		# log prior mass estimate, then scalars
		f = open('output.txt', 'a')
		logw = np.log(frac) + self.log_prior_mass
		line = str(logw)
		for s in self.all_scalars[which, :]:
			line += ' ' + str(s)
		f.write(line + '\n')
		f.close()

		# Forbid another rectangle
		self.forbidden_rectangles = np.vstack([self.forbidden_rectangles, \
										self.all_scalars[which, :]])
		self.refresh_particle(which)

	def refresh_particle(self, which, mcmc_steps=1000):
		"""
		Replace a particle by cloning another and doing MCMC.
		"""
		# Choose a particle to clone
		copy = rng.randint(self.num_particles)
		while copy == which:
			copy = rng.randint(self.num_particles)

		# Clone it
		self.walkers[which] = deepcopy(self.walkers[copy])

		# Do MCMC
		num_accepted = 0
		for i in range(0, mcmc_steps):
			proposal = deepcopy(self.walkers[which])
			logH = proposal.proposal()

			if rng.rand() <= np.exp(logH) and not np.any\
				(Sampler.is_in_rectangle(proposal.scalars, self.forbidden_rectangles)):
				self.walkers[which] = proposal
				num_accepted += 1

	@property
	def rectangle_counts(self):
		"""
		Count how many other walkers are within the rectangle
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
		Is 'scalars' inside the rectangle defined by 'rectangle'?
		Can vectorise on one argument or the other by passing in a 2D array.
		"""
		return np.all(scalars < rectangle, axis=1)


num_particles = 1000
depth = 1000.
sampler = Sampler(num_particles)
sampler.initialise()

for i in range(0, 10):#int(num_particles*depth))
	sampler.do_iteration()

