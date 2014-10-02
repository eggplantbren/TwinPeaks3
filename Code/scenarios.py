from pylab import *

logw = loadtxt('logw.txt')
scalars = loadtxt('scalars.txt')

plot(scalars[:,0], scalars[:,1], 'bo', alpha=0.1)
show()


