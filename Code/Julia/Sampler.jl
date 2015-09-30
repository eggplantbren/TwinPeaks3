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

