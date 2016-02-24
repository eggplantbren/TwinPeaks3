@doc """
An object of this class is a Sampler, obviously.
""" ->
type Sampler
	# The walkers/particles
	num_particles::Int64
	particles::Array{Particle, 1}
	scalars::Array{Float64, 2}

	# Stuff related to the ordering of the scalars
	# by scalar 1 (first column) and then by scalar two (second column)
	indices::Array{Int64, 2}
	ranks::Array{Int64, 2}
end

@doc """
A constructor. Input the number of particles.
""" ->
function Sampler(num_particles::Int64)
	return Sampler(num_particles, Array(Particle, (num_particles, )),
									Array(Float64, (num_particles, 2)),
									Array(Int64, (num_particles, 2)),
									Array(Int64, (num_particles, 2)))
end

@doc """
Generate all the particles from the prior.
""" ->
function initialise!(sampler::Sampler)
	# Generate all the particles from the prior
	# and evaluate their scalars
	for(i in 1:sampler.num_particles)
		sampler.particles[i] = Particle()
		from_prior!(sampler.particles[i])
		sampler.scalars[i, :] = calculate_scalars(sampler.particles[i])
	end

	# Do some sorting.
	sort_scalars!(sampler)
	return nothing
end

@doc """
Sort the scalars. Store various kinds of results
""" ->
function sort_scalars!(sampler::Sampler)
	# Do an argsort
	sampler.indices[:,1] = sortperm(sampler.scalars[:,1])
	sampler.indices[:,2] = sortperm(sampler.scalars[:,2])

	# Use the argsort resuls to compute
	# the ranks of the original array elements
	for(i in 1:sampler.num_particles)
		sampler.ranks[sampler.indices[i, 1], 1] = i
		sampler.ranks[sampler.indices[i, 2], 2] = i
	end
	return nothing
end

@doc """
Calculate the ucc as a function of rank wrt two objective functions.
""" ->
function calculate_uccs(sampler::Sampler)
	uccs = zeros(UInt16, (sampler.num_particles, sampler.num_particles))
	for(i in 1:sampler.num_particles)
		uccs[1:sampler.ranks[i, 1], 1:sampler.ranks[i, 2]] += 1
	end
	# Flip vertically
	uccs = uccs[sampler.num_particles:-1:1, :]
	return uccs
end

@doc """
Calculate the ucc as a function of rank wrt two objective functions.
""" ->
function calculate_uccs2(sampler::Sampler)
	empirical_measure = zeros(Bool, (sampler.num_particles, sampler.num_particles))
	for(i in 1:sampler.num_particles)
		empirical_measure[sampler.ranks[i, 1], sampler.ranks[i, 2]] = true
	end

	uccs = zeros(UInt16, (sampler.num_particles, sampler.num_particles))
	for(i in 1:sampler.num_particles)
		for(j in 1:sampler.num_particles)
			
		end
	end

	return uccs
end

