using PyCall
@pyimport matplotlib.pyplot as plt

include("models/Example.jl")
include("Sampler.jl")

sampler = Sampler(100)
initialise!(sampler)
worst = find_worst(sampler)

plt.scatter(sampler.all_scalars[:,1], sampler.all_scalars[:,2], marker="o",
						s=sampler.uccs)
plt.scatter(sampler.all_scalars[worst,1], sampler.all_scalars[worst,2],
					marker="o", color="r", s=100)
plt.axis("scaled")
plt.axis([0, 1, 0, 1])
plt.show()

