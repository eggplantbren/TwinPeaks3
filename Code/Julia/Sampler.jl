@doc """
An object of this class is a Sampler, obviously.
""" ->
type Sampler
	num_particles::Int64
	particles::Vector{Particle}
	all_scalars::Matrix{Float64}
end

@doc """
A constructor. Input the number of particles.
""" ->
function Sampler(num_particles::Int64)
	return Sampler(num_particles, Array(Particle, (num_particles, )),
									Array(Float64, (num_particles, 2)))
end

@doc """
Generate all the particles from the prior.
""" ->
function initialise!(sampler::Sampler)
	for(i in 1:sampler.num_particles)
		sampler.particles[i] = Particle()
		from_prior!(sampler.particles[i])
		sampler.all_scalars[i, :] = calculate_scalars(sampler.particles[i])
	end
end


