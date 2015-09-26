import numpy as np
import numpy.random as rng

class Walker:
	"""
	Class defining a walker, and the scalar functions
	(log-likelihood or -energy) that you want to increase
	"""
	def __init__(self):
		"""
		Constructor: does nothing.
		"""
		pass

	def from_prior(self):
		"""
		Generate coordinates/parameters from the prior
		"""
		self.x = rng.rand(100)

	def proposal(self):
		"""
		Do something to generate a Metropolis proposal
		"""
		# Choose a probability from log p ~ uniform
		p = 10.**(-3.*rng.rand())

		# Move each coordinate with probability p
		which = rng.rand(len(self.x)) < p

		# Make sure we're at least moving one
		num_moving = which.sum()
		if num_moving == 0:
			which[rng.randint(len(self.x))] = True
			num_moving = 1

		# Move using a heavy-tailed proposal
		self.x[which] += 10.**(1.5 - 6.*rng.rand(num_moving))\
								*rng.randn(num_moving)

		# Mod into prior range
		self.x = np.mod(self.x, 1.)

		# Return zero -- in general return log of
		# a factor that enforces detailed balance
		# wrt the prior.
		return 0.

	@property
	def scalars(self):
		"""
		Calculate the scalar functions of interest
		and return them in the form of a numpy array.
		Current demo problem: the first scalar wants
		the coordinates to cluster in a gaussian in
		the center of the domain, and the second wants
		them to have periodic density. The coordinates
		are IID in this problem so it's easy to numerically
		compute the true partition function.
		"""
		s = np.empty(2)
		s[0] = -0.5*np.sum(self.x**2)
		s[1] = -np.sum(np.sin(4.*np.pi*self.x)**2)
		return s

