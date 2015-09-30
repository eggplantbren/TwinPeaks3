include("Sampler.jl")

srand(0)
sampler = Sampler(100)
initialise!(sampler)

using PyCall
@pyimport matplotlib.pyplot as plt
plt.ion()
plt.hold(true)

while(true)
	scalars = do_iteration!(sampler)
	plt.plot(scalars[1], scalars[2], "b.")
	plt.draw()
end
plt.ioff()
plt.show()

