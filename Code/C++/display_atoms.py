# Display file for Gravity
from pylab import *

output = atleast_2d(loadtxt('sample.txt'))

ion()
hold(False)
for i in xrange(0, output.shape[0]):
	model = output[i, 1:]
	x, y = model[0:50], model[100:150]
	plot(x, y, 'k.')
	axis('scaled')
	axis([0, 1, 0, 1])
	title(i+1)
	draw()

ioff()
show()

