include("Sampler.jl")

srand(0)
sampler = Sampler(1000)
initialise!(sampler)

using PyCall
@pyimport matplotlib.pyplot as plt
plt.ion()
plt.hold(true)

for(i in 1:10000000)
	scalars = do_iteration!(sampler)

	plt.plot(scalars[1], scalars[2], "b.", markersize=1)
	if(rem(i, 100) == 0)
		plt.draw()
	end
end
plt.ioff()
plt.show()
