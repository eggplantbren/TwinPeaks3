# Display file for Gravity
from pylab import *

output = atleast_2d(loadtxt('sample.txt'))

ion()
hold(False)
for i in xrange(0, output.shape[0]):
	model = output[i, 1:]
	x, y = model[0:1000], model[1000:2000]
	plot(x, y, 'k.')
	axis('scaled')
	title(i+1)
	draw()

ioff()
show()

