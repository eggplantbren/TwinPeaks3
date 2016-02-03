# Display file for Potts
from pylab import *

output = atleast_2d(loadtxt('sample.txt'))

ion()
hold(False)
for i in range(0, output.shape[0]):
	x = output[i, 3:-1].reshape((100, 100))
	imshow(x, interpolation='nearest')
	title(str(i+1) + '/' + str(output.shape[0]))
	print(output[i, -1])
	draw()
ioff()
show()

