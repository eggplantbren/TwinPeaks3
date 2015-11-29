# Display file for Gravity
from pylab import *

output = atleast_2d(loadtxt('sample.txt'))
num_particles = (output.shape[1] - 3)//6
print(num_particles)

figure(figsize=(13, 6))
ion()
hold(False)
for i in xrange(0, output.shape[0]):
	model = output[i, 3:-1]
	x, y = model[0:num_particles], model[num_particles:2*num_particles]
	vx, vy = model[3*num_particles:4*num_particles], model[4*num_particles:5*num_particles]
	subplot(1,2,1)
	plot(x, y, 'k.')
	axis([-10, 10, -10, 10])
	subplot(1,2,2)
	plot(vx, vy, 'k.')
	axis([-10, 10, -10, 10])
	title(i+1)
	draw()

ioff()
show()

