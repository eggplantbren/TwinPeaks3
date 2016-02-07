"""
Experimenting with some animations to provide intuition about
the TwinPeaks algorithm.
\pi(X1, X2) ~ Uniform([0, 1]^2)
Constraint on X1X2 < const
const decreases.
"""

import numpy as np
import matplotlib.pyplot as plt

def deriv(state):
	return (state - state*np.log(state))/np.log(state)

# Let x=X1, y=X2|X1
x = np.linspace(1E-6, 1, 1001)

# Timestep
t_final, dt = 10., 0.001
steps = int(t_final/dt)

# Initial condition for threshold
C = 0.99

# Set up plotting
plt.figure(figsize=(8, 8))
plt.ion()
plt.hold(False)

for i in range(0, steps):
	f1 = deriv(C)
	f2 = deriv(C + 0.5*dt*f1)
	f3 = deriv(C + 0.5*dt*f2)
	f4 = deriv(C + dt*f3)
	C += dt/6*(f1 + 2*f2 + 2*f3 + f4)

	y = C/x
	plt.plot(x, y, 'k-', linewidth=2)
	plt.xlabel('$x$', fontsize=18)
	plt.ylabel('$y$', fontsize=18)
	plt.axis([0, 1, 0, 1])
	plt.title('$X = {X:.4f}$'.format(X=C-C*np.log(C)), fontsize=16)

	yy = y.copy()
	yy[y > 1.] = 1.
	plt.fill_between(x, yy, color='cyan', alpha=0.3)
	plt.draw()

# Finalise plotting
plt.ioff()
plt.show()
