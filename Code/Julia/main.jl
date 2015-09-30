include("Sampler.jl")

sampler = Sampler(100)
initialise!(sampler)

using PyCall
@pyimport matplotlib.pyplot as plt

counts = rectangle_counts(sampler)
plt.hist(counts, 100)
plt.show()

do_iteration!(sampler)

