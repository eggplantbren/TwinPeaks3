@doc """
An object of this class is a Sampler, obviously.
""" ->
type Sampler
	num_particles::Int64
	particles::Vector{Particle}
	all_scalars::Matrix{Float64}
	uccs::Vector{Int64}
	tiebreakers::Vector{Float64}
end

@doc """
A constructor. Input the number of particles.
""" ->
function Sampler(num_particles::Int64)
	return Sampler(num_particles, Array(Particle, (num_particles, )),
									Array(Float64, (num_particles, 2)),
									Array(Int64, (num_particles, )),
									Array(Float64, (num_particles, )))
end

@doc """
Generate all the particles from the prior.
""" ->
function initialise!(sampler::Sampler)
	for(i in 1:sampler.num_particles)
		sampler.particles[i] = Particle()
		from_prior!(sampler.particles[i])
		sampler.all_scalars[i, :] = calculate_scalars(sampler.particles[i])
		sampler.tiebreakers[i] = rand()
	end
	calculate_uccs!(sampler::Sampler)
end

@doc """
Calculate the upper corner count of each particle.
""" ->
function calculate_uccs!(sampler::Sampler)
	for(i in 1:sampler.num_particles)
		sampler.uccs[i] = 0
		for(j in 1:sampler.num_particles)
			if(all(sampler.all_scalars[j, :] .> sampler.all_scalars[i, :]))
				sampler.uccs[i] += 1
			end
		end
	end
end

@doc """
Find the particle with the highest (ucc, tiebreaker) combo.
""" ->
function find_worst(sampler::Sampler)
	which = 1
	for(i in 2:sampler.num_particles)
		worse = false
		if(sampler.uccs[i] > sampler.uccs[which])
			worse = true
		end
		if(sampler.uccs[i] == sampler.uccs[which])
			if(sampler.tiebreakers[i] > sampler.tiebreakers[which])
				worse = true
			end
		end

		if(worse)
			which = i
		end
	end
	return which
end

@doc """
Do one iteration of Nested Sampling.
""" ->
function do_iteration!(sampler::Sampler)
	worst = find_worst(sampler)
	println(sampler.uccs[worst], " ", sampler.tiebreakers[worst])
end
