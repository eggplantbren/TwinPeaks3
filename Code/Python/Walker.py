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
		"""
		pass

