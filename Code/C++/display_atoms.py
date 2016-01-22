# Display file for Atoms
from pylab import *
import os
from mpl_toolkits.mplot3d import Axes3D

save_plots = True
if save_plots:
	os.system('rm Frames/*.png')

output = atleast_2d(loadtxt('sample.txt'))
num_atoms = (output.shape[1] - 3)//3

fig = figure(figsize=(12, 9))
ion()
ax = fig.add_subplot(111, projection='3d')
hold(False)
k = 0
for i in range(0, output.shape[0]):
	model = output[i, 3:]
	x, y, z = model[0:num_atoms], model[num_atoms:2*num_atoms],\
						model[2*num_atoms:3*num_atoms]
	ax.scatter(x, y, z, marker='o', s=20, c='goldenrod', alpha=0.6)
	axis('scaled')
	axis([0, 100, 0, 100])
	title(i+1)
	xlabel('$x$', fontsize=16)
	ylabel('$y$', fontsize=16)
	draw()

	if save_plots:
		savefig('Frames/' + '%0.4d'%(k+1) + '.png', bbox_inches='tight')
		print('Frames/' + '%0.4d'%(k+1) + '.png')
	k += 1

# Now rotate the camera
for ii in linspace(0, 360, 101):
	ax.scatter(x, y, z, marker='o', s=20, c='goldenrod', alpha=0.6)
	ax.view_init(elev=30., azim=ii)
	draw()
	if save_plots:
		savefig('Frames/' + '%0.4d'%(k+1) + '.png', bbox_inches='tight')
		print('Frames/' + '%0.4d'%(k+1) + '.png')
	k += 1

ioff()
show()

