from pylab import *

output = atleast_2d(loadtxt('sample.txt'))

figure(figsize=(13, 6))
ion()
hold(False)
for i in range(0, output.shape[0]):
	model = output[i, 3:]
	hist(model, 300, alpha=0.2)
	axis([0., 1., 0., 10.])
	title(i+1)
	draw()

ioff()
show()

