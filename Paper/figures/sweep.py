"""
Experimenting with some animations to provide intuition about
the TwinPeaks algorithm.
\pi(X1, X2) ~ Uniform([0, 1]^2)
Constraint on X1X2 < const
const decreases.
"""

import numpy as np
import matplotlib.pyplot as plt

# Let x=X1, y=X2|X1

dt = 0.01

for i in range(0, steps):

