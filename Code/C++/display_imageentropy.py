# Display file for Potts
from pylab import *

output = atleast_2d(loadtxt('sample.txt'))

ion()
hold(False)
for i in xrange(0, output.shape[0]):
	x = output[i, 1:].reshape((100, 100))
	imshow(x, interpolation='nearest', cmap='gray')
	title(str(i+1) + '/' + str(output.shape[0]))
	draw()
ioff()
show()

