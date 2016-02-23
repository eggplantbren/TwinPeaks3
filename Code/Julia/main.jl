# Setup plotting and include other files
using PyCall
@pyimport matplotlib.pyplot as plt
@pyimport colormaps
include("Models/Example.jl")
include("Sampler.jl")

# Create and initialise a Sampler
sampler = Sampler(1000)
initialise!(sampler)

# Plot scalars and ranks
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.scatter(sampler.scalars[:, 1], sampler.scalars[:, 2], marker="o", alpha=0.1)

plt.subplot(1, 2, 2)
plt.scatter(sampler.ranks[:,1], sampler.ranks[:,2], marker="o", alpha=0.1)
plt.axis([0, 1+sampler.num_particles+1, 0, 1+sampler.num_particles])
plt.show()

uccs = calculate_uccs(sampler)
plt.imshow(uccs, interpolation="nearest", cmap="Blues")
plt.show()

