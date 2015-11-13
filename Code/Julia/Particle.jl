@doc """
Particle class
""" ->
type Particle
	params::Array{Float64, 1}
end

@doc """
Constructor taking no arguments
""" ->
function Particle()
	return Particle()
end

