using PyCall
@pyimport matplotlib.pyplot as plt

include("models/Example.jl")
include("Sampler.jl")

sampler = Sampler(1000)
initialise!(sampler)
indices = sortperm(sampler.uccs)
reverse!(indices)

plt.ion()
plt.hold(true)
plt.scatter(sampler.all_scalars[:,1], sampler.all_scalars[:,2], marker="o",
						s=0.5*sampler.uccs, alpha=0.1)
plt.axis("scaled")
plt.axis([0, 1, 0, 1])

med = median(sampler.uccs)
for(i in 1:div(sampler.num_particles, 2))
	ii = indices[i]

	plt.fill_between([0.0, sampler.all_scalars[ii, 1]],
						sampler.all_scalars[ii, 2], color=[0.5, 0.5, 0.5])
	plt.title(string(i, " rectangles."))
	plt.draw()
end

plt.ioff()
plt.show()

