@doc """
An object of this class is a Sampler, obviously.
""" ->
type Sampler
	num_particles::Int64
	particles::Vector{Particle}
	scalars::Array{Float64, 2}
	scalars_sorted::Array{Float64, 2}
end

@doc """
A constructor. Input the number of particles.
""" ->
function Sampler(num_particles::Int64)
	return Sampler(num_particles, Array(Particle, (num_particles, )),
									Array(Float64, (num_particles, 2)),
									Array(Float64, (num_particles, 2)))
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

	# Sort the scalars
	sampler.scalars_sorted[:,1] = sort(sampler.scalars[:,1])
	sampler.scalars_sorted[:,2] = sort(sampler.scalars[:,2])
end

@doc """
Calculate the ucc at a certain position (defined by a pair of ranks)
""" ->
function calculate_ucc(sampler::Sampler, ranks::Tuple{Int64, Int64})
	@assert length(ranks) == 2

	# Scalars corresponding to the pair of ranks
	s = (sampler.scalars_sorted[ranks[1], 1],
					sampler.scalars_sorted[ranks[2], 2])

	
end

