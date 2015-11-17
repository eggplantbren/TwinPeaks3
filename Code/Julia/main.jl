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
rectangles = Array(Float64, (sampler.num_particles, 2))

for(i in 1:div(sampler.num_particles, 2))
	ii = indices[i]
	rectangles[i, :] = sampler.all_scalars[ii, :]
end
rectangles = rectangles[1:div(sampler.num_particles, 2), :]

for(i in 1:size(rectangles)[1])
	plt.fill_between([0.0, rectangles[i, 1]],
						rectangles[i, 2], color=[0.5, 0.5, 0.5])
	plt.title(string(i, " rectangles."))
	if(rem(i, 10) == 0)
		plt.draw()
	end
end

# Estimate the fraction of prior mass we just excluded
reps = 100000
num_ok = 0
for(i in 1:reps)
	particle = Particle()
	from_prior!(particle)
	scalars = calculate_scalars(particle)
	is_okay = true
	for(j in 1:size(rectangles)[1])
		if(all(scalars .< vec(rectangles[j, :])))
			is_okay = false
		end
	end
	if(is_okay)
		num_ok += 1
		plt.plot(calculate_scalars(particle)[1], calculate_scalars(particle)[2], "r.", markersize=3)
	end
	if(rem(i, 100) == 0)
		plt.draw()
	end
	println(i, " ", num_ok)
end

plt.ioff()
plt.show()
