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

walker = Walker()
from_prior!(walker)

