# Setup plotting and include other files
using PyCall
@pyimport matplotlib.pyplot as plt
include("models/Example.jl")
include("Sampler.jl")

# Use an odd number of particles
num_particles = 101
@assert mod(num_particles, 2) != 0

# Create and initialise a Sampler
sampler = Sampler(num_particles)
initialise!(sampler)

# argsort the uccs
indices = sortperm(sampler.uccs)
reverse!(indices)

# Find the ucc threshold
threshold_index = 1 + div(num_particles, 2)
while((threshold_index != num_particles)
		&& (sampler.uccs[threshold_index+1] == sampler.uccs[threshold_index]))
	threshold_index += 1
end

println(threshold_index)
#threshold = sampler.uccs[indices[]]

plt.scatter(sampler.all_scalars[:,1], sampler.all_scalars[:,2], marker="o",
						s=0.5*sampler.uccs, alpha=0.1)
plt.axis("scaled")
plt.axis([0, 1, 0, 1])
plt.show()

