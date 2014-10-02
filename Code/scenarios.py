from pylab import *
import pylab

logw = loadtxt('logw.txt')
scalars = loadtxt('scalars.txt')

plot(scalars[:,0], scalars[:,1], 'bo', alpha=0.1)
show()

# Truncate so we only use completed runs
steps = 200
runs = logw.size//steps

scalars = scalars[0:runs*steps, :]

# Generate logX values that have a valid ordering.
# Strictly decreasing within each run, and any point
# with BOTH scalars higher than another point must have
# lower logX.

logX = zeros(scalars.shape[0])

# Make an arbitrary assignment for run 1
for i in xrange(0, steps):
  logX[i] = -(i+1)

# Go through other runs
for i in xrange(steps, runs*steps):
  lower_limit, upper_limit = -1E300, 0.

  if i % steps != 0:
    upper_limit = logX[i-1]

  # Go through all assigned values in previous runs, and find
  # all distributions known to be exterior to the current one
  # we're trying to assign. If any have an assignment below
  # upper_limit, lower the upper_limit.
  for j in xrange(0, (i//steps)*steps):
    # If exterior
    if pylab.all(scalars[j, :] < scalars[i, :]) and logX[j] < upper_limit:
      upper_limit = logX[j]

  # Do the same for the lower limit. Look for any distributions
  # known to be interior to the current point. Raise lower_limit
  # if necessary.
  for j in xrange(0, (i//steps)*steps):
    # If exterior
    if pylab.all(scalars[j, :] > scalars[i, :]) and logX[j] > lower_limit:
      lower_limit = logX[j]

  if lower_limit == -1E300:
    logX[i] = logX[0:i].min() - 1
  else:
    logX[i] = 0.5*(lower_limit + upper_limit)

plot(logX, scalars[:,0], 'bo')
show()

