# Setup plotting and include other files
using PyCall
@pyimport matplotlib.pyplot as plt
@pyimport colormaps
include("Models/Example.jl")
include("Sampler.jl")

# Create and initialise a Sampler
sampler = Sampler(1000)
initialise!(sampler)

# Plot scalars
plt.scatter(sampler.scalars[:, 1], sampler.scalars[:, 2],
					marker="o", alpha=0.5)
plt.axis("scaled")
plt.axis([-0.01, 1.01, -0.01, 1.01])
plt.show()

