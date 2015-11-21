# Setup plotting and include other files
using PyCall
@pyimport matplotlib.pyplot as plt
@pyimport colormaps
include("models/Example.jl")
include("Sampler.jl")

# Use an odd number of particles
num_particles = 101
@assert mod(num_particles, 2) != 0

# Create and initialise a Sampler
sampler = Sampler(num_particles)
initialise!(sampler)

# argsort the uccs from highest to lowest
indices = sortperm(sampler.uccs)
reverse!(indices)

# Find the ucc threshold
threshold_index = 1 + div(num_particles, 2)
while((threshold_index != num_particles)
		&& (sampler.uccs[indices[threshold_index+1]] == sampler.uccs[indices[threshold_index]]))
	threshold_index += 1
end
threshold = sampler.uccs[indices[threshold_index]]

# Which particles will die?
dying = sampler.uccs .>= threshold

# Plot scalars
plt.scatter(sampler.all_scalars[dying,1], sampler.all_scalars[dying,2],
					marker="o", alpha=0.5)
plt.axis("scaled")
plt.axis([-0.01, 1.01, -0.01, 1.01])
plt.show()

