from pylab import *

x = linspace(0., 1., 2001)
p = exp(-10.*(x - 0.5)**2 - 1.*sin(2.*pi*x/0.5)**2)
Z = trapz(p, x=x)
H = trapz(p/Z*log(p/Z), x=x)
print(200*log(Z))
print(200*H)
plot(x, p/Z)
show()

