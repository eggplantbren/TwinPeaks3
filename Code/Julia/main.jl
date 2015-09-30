include("Sampler.jl")

srand(0)
sampler = Sampler(100)
initialise!(sampler)

using PyCall
@pyimport matplotlib.pyplot as plt
plt.ion()
plt.hold(true)

for(i in 1:5000)
	scalars = do_iteration!(sampler)

	plt.plot(scalars[1], scalars[2], "b.")
	if(rem(i, 10) == 0)
		plt.draw()
	end
end
plt.ioff()
plt.show()

