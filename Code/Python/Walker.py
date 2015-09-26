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
		self.x = rng.rand(10)

	def proposal(self):
		"""
		Do something to generate a Metropolis proposal
		"""
		pass

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
		s = np.array(2)
		s[0] = -0.5*np.sum(self.x**2)
		s[1] = -np.sum(np.sin(4.*np.pi*x)**2)
		return s

