include("Walker.jl")

using PyCall
@pyimport matplotlib.pyplot as plt

# An object of this class is a sampler.
type Sampler
	# Number of walkers
	num_walkers::Int64

	# Walkers
	walkers::Array{Walker, 1}

	# Fraction of prior mass still in play
	log_prior_mass::Float64

	# Iteration
	iteration::Int64

	# Forbidden rectangles
	forbidden_rectangles::Array{Float64, 2}
end

# A constructor that just takes the number of walkers
function Sampler(num_walkers::Int64)
	return Sampler(num_walkers, Array(Walker, num_walkers), 0.0, 0,
				Array(Float64, (num_scalars, 0)))
end

# Generate all the walkers from the prior
function initialise!(sampler::Sampler)
	for(i in 1:sampler.num_walkers)
		sampler.walkers[i] = Walker()
		from_prior!(sampler.walkers[i])
	end

	# Open and close output file to clear it
	f = open("output.txt", "w")
	close(f)

	return nothing
end

# Remove redundant rectangles (only works well if run each iteration)
function prune_rectangles!(sampler::Sampler)
	n = size(sampler.forbidden_rectangles)[2]
	keep = fill(true, n)
	for(i in 1:(n-1))
		if(is_in_rectangle(sampler.forbidden_rectangles[:,i],
							sampler.forbidden_rectangles[:,n]))
			keep[i] = false
		end
	end
	sampler.forbidden_rectangles = sampler.forbidden_rectangles[:, keep]
end

# Do an NS iteration
function do_iteration!(sampler::Sampler)
	counts = rectangle_counts(sampler)

	# Choose one with minimum count to discard
	temp = find(counts .== 0)
	which = temp[rand(1:length(temp))]

	# Retain its scalars
	scalars = sampler.walkers[which].scalars

	# Estimate the fraction of *remaining* prior mass being eliminated
	frac = (1 + counts[which])/sampler.num_walkers

	# Write log(prior mass) and scalars
	f = open("output.txt", "a")
	print(f, log(frac) + sampler.log_prior_mass, " ")
	for(s in scalars)
		print(f, s, " ")
	end
	print(f, "\n")
	close(f)

	# Reduce prior mass remaining
	sampler.log_prior_mass += log(1.0 - frac)

	# Forbid another rectangle
	forbid_rectangle!(sampler, scalars)
	prune_rectangles!(sampler)

	# Refresh discarded walker
	refresh_walker!(sampler, which)

	return scalars
end

# Add an extra rectangle to the list of forbidden ones
function forbid_rectangle!(sampler::Sampler, scalars::Array{Float64, 1})
	sampler.forbidden_rectangles = hcat(sampler.forbidden_rectangles, scalars)
	return nothing
end

# Is 'scalars' in the allowed region?
function is_okay(sampler::Sampler, scalars::Array{Float64, 1})
	for(i in 1:size(sampler.forbidden_rectangles)[2])
		if(is_in_rectangle(scalars, sampler.forbidden_rectangles[:,i]))
			return false
		end
	end
	return true
end

# Do MCMC to replace one of the walkers
function refresh_walker!(sampler::Sampler, which::Int64,
									mcmc_steps::Int64=1000)
	# Choose a particle to clone
	copy = rand(1:sampler.num_walkers)
	while copy == which
		copy = rand(1:sampler.num_walkers)
	end

	# Clone it
	sampler.walkers[which] = deepcopy(sampler.walkers[copy])

	# Do MCMC
	num_accepted = 0
	for(i in 1:mcmc_steps)
		proposal = deepcopy(sampler.walkers[which])
		logH = proposal!(proposal)

		if((rand() <= exp(logH)) & is_okay(sampler, proposal.scalars))
			sampler.walkers[which] = proposal
			num_accepted += 1
		end
	end

	println("Iteration ", sampler.iteration+1, ", accepted ", num_accepted, "/",
						mcmc_steps)

	sampler.iteration += 1
end

# Count how many other walkers are within the rectangle of each walker.
function rectangle_counts(sampler::Sampler)
	counts = Array(Int64, sampler.num_walkers)
	for(i in 1:sampler.num_walkers)
		counts[i] = 0
		for(j in 1:sampler.num_walkers)
			counts[i] += is_in_rectangle(sampler.walkers[j].scalars,
											sampler.walkers[i].scalars)
		end
	end
	return counts
end

# Is 'scalars' inside the rectangle defined by 'rectangle'?
function is_in_rectangle(scalars::Array{Float64, 1},
							rectangle::Array{Float64, 1})
	for(i in 1:length(scalars))
		if(scalars[i] >= rectangle[i])
			return false
		end
	end
	return true
end
