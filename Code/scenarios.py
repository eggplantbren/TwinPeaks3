from pylab import *

# Generate from a negative-exponential
def generate_exponential(lower, upper):
  u_max = 1. - exp(-(upper - lower))
  if isinf(lower):
    u = rand()
  else:
    u = rand()*u_max
  return upper + log(1. - u)

# Find which elements of 'others' are immediately adjacent
# to 'logl', assign a logx value between the corresponding limits
def sandwich(scalars, other_scalars, other_logx, upper2):
  above = empty(other_scalars.shape[0]).astype('bool')
  below = empty(other_scalars.shape[0]).astype('bool')
  for i in range(0, above.size):
    above[i] = logical_and(other_scalars[i,0] > scalars[0], other_scalars[i, 1] > scalars[1])
    below[i] = logical_and(other_scalars[i,0] < scalars[0], other_scalars[i, 1] < scalars[1])

  try:
    upper = (other_logx[below]).min()
  except:
    upper = 0.

  try:
    lower = (other_logx[above]).max()
  except:
    lower = -Inf

  upper = array([upper, upper2]).min()

  return generate_exponential(lower, upper)

logw = loadtxt('logw.txt')
scalars = loadtxt('scalars.txt')

plot(scalars[:,0], scalars[:,1], 'bo', alpha=0.1)
show()

# Calculate size of inputs
reps = int(sum(logw == logw.min()))
steps = len(logw)//reps
logX = empty(logw.size)

# First rep
logX[0] = log(rand())
for i in xrange(1, steps):
  logX[i] = logX[i-1] + log(rand())

# The rest
for i in xrange(steps, logw.size):
  if i % steps == 0:
    upper2 = 0.
  else:
    upper2 = logX[i-1]
  logX[i] = sandwich(scalars[i, :], scalars[0:i, :], logX[0:i], upper2)

plot(logX, 'b')
savetxt('logw2.txt', logX)
show()

