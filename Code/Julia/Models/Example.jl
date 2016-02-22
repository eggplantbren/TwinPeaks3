include("../Utils.jl")

@doc """
An object of this class represents a point in parameter space.
There are functions defined to evaluate the log likelihood and
move around.
""" ->
type Particle
	params::Vector{Float64}
end

@doc """
A constructor taking no parameters.
""" ->
function Particle()
	return Particle(Array(Float64, (2, )))
end

@doc """
Generate params from the prior
""" ->
function from_prior!(particle::Particle)
	particle.params[1] = randn()
	particle.params[2] = particle.params[1] + randn()
	return nothing
end

@doc """
Do a metropolis proposal. Return log(hastings factor for prior sampling)
""" ->
function perturb!(particle::Particle)
	log_H = 0.0

	which = rand(1:2)
	if(which == 1)
		log_H -= -0.5*particle.params[1]^2
		particle.params[1] += randh()
		log_H += -0.5*particle.params[1]^2
	else
		particle.params[2] -= particle.params[1]
		log_H -= -0.5*particle.params[2]^2
		particle.params[2] += randh()
		log_H += -0.5*particle.params[2]^2
		particle.params[2] += particle.params[1]
	end
	return log_H
end

@doc """
Evaluate the two objective functions
""" ->
function calculate_scalars(particle::Particle)
	s = Array(Float64, (2, ))
	s[1] = particle.params[1]
	s[2] = particle.params[2]
	return s
end

@doc """
Convert to string, for output to sample.txt
"""
import Base.string
function string(particle::Particle)
	return join([string(signif(x, 6), " ") for(x in particle.params)])
end

