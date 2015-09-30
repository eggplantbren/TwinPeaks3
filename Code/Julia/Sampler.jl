include("Walker.jl")

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
	return nothing
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

	# Reduce prior mass remaining
	sampler.log_prior_mass += log(1.0 - frac)

	# Forbid another rectangle
	forbid_rectangle!(sampler, scalars)

	# Refresh discarded walker
	refresh_walker!(sampler, which)

	return scalars
end

# Add an extra rectangle to the list of forbidden ones
function forbid_rectangle!(sampler::Sampler, scalars::Array{Float64, 1})
	return nothing
end

# Do MCMC to replace one of the walkers
function refresh_walker!(sampler::Sampler, which::Int64)
	return nothing
end


# Count how many other walkers are within the rectangle of each walker.
function rectangle_counts(sampler::Sampler)
	counts = Array(Int64, sampler.num_walkers)
	for(i in 1:sampler.num_walkers)
		counts[i] = 0
		for(j in 1:sampler.num_walkers)
			counts[i] += is_in_rectangle(sampler.walkers[i].scalars,
											sampler.walkers[j].scalars)
		end
	end
	return counts
end

# Is 'scalars' inside the rectangle defined by 'rectangle'?
function is_in_rectangle(scalars::Array{Float64, 1},
							rectangle::Array{Float64, 1})
	return all(scalars .< rectangle)
end

