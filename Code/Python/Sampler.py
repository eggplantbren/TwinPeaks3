import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
from copy import deepcopy

from Particle import Particle

class Sampler:
	"""
	An object of this class is a sampler.
	"""

	def __init__(self, num_particles):
		"""
		Constructor: pass in the number of particles
		"""
		self.num_particles = num_particles
		self.particles = [Particle() for i in range(0, num_particles)]

		# Open and close output file to clear it
		f = open('output.txt', 'w')
		f.close()

		# Rectangles that have been forbidden so far
		self.forbidden_rectangles = None

		# Fraction of prior mass still in play
		self.log_prior_mass = 0.

	def initialise(self):
		"""
		Generate all the particles from the prior
		"""
		self.all_scalars = []	# Values of scalars for all particles
		for particle in self.particles:
			particle.from_prior()
			self.all_scalars.append(particle.scalars)
		self.all_scalars = np.array(self.all_scalars)
		self.forbidden_rectangles = np.empty((0, self.all_scalars.shape[1]))
		self.iteration = 0
		return

	def do_iteration(self):
		"""
		Do an NS iteration
		"""
		counts = self.corner_counts

		# Choose one with minimum count to discard
		temp = np.nonzero(counts == counts.min())[0]
		which = temp[rng.randint(len(temp))]

		# Need to figure out what happens in this case
		if counts.min() != 0:
			print("counts.min() != 0")

		# Estimate the fraction of *remaining* prior mass being eliminated
		frac = float(1 + counts[which])/self.num_particles

		# Copy out the scalars for returning purposes
		keep = self.all_scalars[which, :].copy()

		# Write discarded particle info to disk
		# log prior mass estimate, then scalars
		f = open('output.txt', 'a')
		logw = np.log(frac) + self.log_prior_mass
		line = str(logw)
		for s in self.all_scalars[which, :]:
			line += ' ' + str(s)
		f.write(line + '\n')
		f.close()

		# Reduce prior mass remaining
		self.log_prior_mass += np.log(1. - frac)

		# Forbid another rectangle
		self.forbidden_rectangles = np.vstack([self.forbidden_rectangles, \
										self.all_scalars[which, :]])
		self.refresh_particle(which)

		return keep

	def refresh_particle(self, which, mcmc_steps=1000):
		"""
		Replace a particle by cloning another and doing MCMC.
		"""
		# Choose a particle to clone
		copy = rng.randint(self.num_particles)
		while copy == which:
			copy = rng.randint(self.num_particles)

		# Clone it
		self.particles[which] = deepcopy(self.particles[copy])

		# Do MCMC
		num_accepted = 0
		for i in range(0, mcmc_steps):
			proposal = deepcopy(self.particles[which])
			logH = proposal.proposal()
			proposal_scalars = proposal.scalars

			if rng.rand() <= np.exp(logH) and not np.any\
				(Sampler.is_in_rectangle\
					(proposal_scalars, self.forbidden_rectangles)):
				self.particles[which] = proposal
				self.all_scalars[which, :] = proposal_scalars
				num_accepted += 1

		print('Iteration {i}. Accepted {a}/{b}.'.format(i=(self.iteration+1), a=num_accepted, b=mcmc_steps))
		self.iteration += 1

	@property
	def corner_counts(self):
		"""
		Count how many other particles are within the rectangle
		of each particle.
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

if __name__ == '__main__':
	num_particles = 100
	depth = 50.			# How far to go
	sampler = Sampler(num_particles)
	sampler.initialise()

	plt.ion()
	plt.hold(True)

	for i in range(0, int(num_particles*depth)):
		keep = sampler.do_iteration()
		plt.plot(keep[0], keep[1], 'b.')
		plt.xlabel('Scalar 1')
		plt.ylabel('Scalar 2')
		plt.title('Discarded Points')
		plt.draw()

	plt.ioff()
	plt.show()

