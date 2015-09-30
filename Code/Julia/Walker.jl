include("Utils.jl")

# Walker class. Defines what problem we're solving.
# An object of this class is a point in parameter space. 
type Walker
	# Number of dimensions
	N::Int64

	# The parameters
	params::Array{Float64, 1}

	# The scalars
	scalars::Array{Float64, 1}
end

# Constructor that takes no input
function Walker()
	# Set the number of dimensions, the parameters, and the scalars
	Walker(100, zeros(100), zeros(2))
end

# Initialise a walker from the prior
function from_prior!(walker::Walker)
	walker.params = rand(walker.N)
	calculate_scalars!(walker)
end

# Calculate the scalars
function calculate_scalars!(walker::Walker)
	walker.scalars[1] = -sum((walker.params - 0.5).^2)
	walker.scalars[2] = -sum(4.*pi*walker.params.^2)
end

# Move a Walker according to the prior
function proposal!(walker::Walker)
	# Do something to generate a Metropolis proposal
	# Choose a probability from log p ~ uniform
	p = 10.0^(-3.0*rand())

	# Move each coordinate with probability p
	which = rand(walker.N) .< p

	# Make sure we're at least moving one
	num_moving = sum(which)
	if num_moving == 0
		which[rand(1:walker.N)] = true
		num_moving = 1
	end

	# Move using a heavy-tailed proposal
	walker.params[which] += randh(num_moving)

	# Mod into prior range
	walker.params = mod(walker.params, 1.0)

	# Calculate scalars
	calculate_scalars!(walker)

	# Return zero.
	# In general return log of a factor that enforces detailed balance
	# wrt the prior.
	return 0.0
end


