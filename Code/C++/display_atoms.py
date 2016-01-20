# Display file for Gravity
from pylab import *

output = atleast_2d(loadtxt('sample.txt'))
num_atoms = (output.shape[1] - 3)//3

ion()
hold(False)
for i in range(0, output.shape[0]):
	model = output[i, 3:]
	x, y = model[0:num_atoms], model[2*num_atoms:3*num_atoms]
	plot(x, y, 'k.')
	axis('scaled')
	axis([0, 10, 0, 10])
	title(i+1)
	draw()

ioff()
show()

