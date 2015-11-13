using PyCall
@pyimport matplotlib.pyplot as plt

include("models/Example.jl")
include("Sampler.jl")

sampler = Sampler(100)
initialise!(sampler)

plt.scatter(sampler.all_scalars[:,1], sampler.all_scalars[:,2], marker="o",
						s=sampler.uccs)
plt.axis("scaled")
plt.axis([0, 1, 0, 1])
plt.show()

