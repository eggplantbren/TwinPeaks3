#include <iostream>
#include <cassert>
#include "Utils.h"

namespace TwinPeaks
{

template<class MyModel>
Sampler<MyModel>::Sampler(unsigned int num_particles, unsigned int mcmc_steps,
                            const std::vector<RNG>& rngs)
:num_particles(num_particles)
,mcmc_steps(mcmc_steps)
,rngs(rngs)
,particles(num_particles)
,scalars(num_scalars, std::vector<ScalarType>(num_particles))
,indices(num_scalars, std::vector<size_t>(num_particles))
,ranks(num_scalars, std::vector<size_t>(num_particles))
,particle_uccs(num_particles)
{
    assert(num_particles > 0 && mcmc_steps > 0);
}

template<class MyModel>
void Sampler<MyModel>::initialise()
{
    std::cout<<"# Initialising TwinPeaks sampler."<<std::endl;
    std::cout<<"# Generating "<<num_particles<<" particles from the prior...";
    std::cout<<std::flush;
    for(unsigned int i=0; i<num_particles; ++i)
    {
        particles[i].from_prior(rngs[0]);
        const std::vector<double>& s = particles[i].get_scalars();
        scalars[i].clear();
        for(size_t j=0; j<s.size(); j++)
            scalars[i].push_back(ScalarType(s[j]));
        for(size_t j=0; j<scalars[i].size(); j++)
            scalars[i][j].from_prior(rngs[0]);
    }
    std::cout<<"done."<<std::endl<<std::endl;
}

template<class MyModel>
void Sampler<MyModel>::do_iteration()
{
    // Argsort by the two scalars
    for(size_t i=0; i<scalars.size(); ++i)
    {
        indices[i] = argsort(scalars[i]);
        ranks[i] = compute_ranks(indices[i]);
    }
}

/*
template<class MyModel>
void Sampler<MyModel>::prune_rectangles(const std::vector<ScalarType>& latest)
{
	// Loop over all rectangles
	for(auto it=rects.begin(); it != rects.end(); ++it)
	{
		if(it->compare(latest) == 1)
			it = rects.erase(it);
	}
}

template<class MyModel>
double Sampler<MyModel>::do_iteration()
{
	// Extract and sort the two scalars
	std::vector<ScalarType> s1(num_particles);
	std::vector<ScalarType> s2(num_particles);
	for(int i=0; i<num_particles; i++)
	{
		s1[i] = scalars[i][0];
		s2[i] = scalars[i][1];
	}
	// Calculate ranks
	std::vector<size_t> r1 = ranks(s1);
	std::vector<size_t> r2 = ranks(s2);
	// Sort
	std::sort(s1.begin(), s1.end());
	std::sort(s2.begin(), s2.end());

	// Calculate empirical measure of ranks
	std::vector< std::vector<bool> > empirical_measure(num_particles,
										std::vector<bool>(num_particles, 0));
	for(int i=0; i<num_particles; i++)
		empirical_measure[num_particles - r2[i] - 1][r1[i]] = 1;

	// Integrate empirical measure to get (inclusive) upper corner count
	std::vector< std::vector<unsigned short> > ucc(num_particles,
										std::vector<unsigned short>(num_particles, 0));
	int up, right, up_and_right;
	for(int i=0; i<num_particles; i++)
	{
		for(int j=(num_particles-1); j>=0; j--)
		{
			up = 0;
			right = 0;
			up_and_right = 0;
			if(i != 0)
				up = ucc[i-1][j];
			if(j != (num_particles-1))
				right = ucc[i][j+1];
			if((i != 0) && (j != (num_particles-1)))
				up_and_right = ucc[i-1][j+1];
			ucc[i][j] = empirical_measure[i][j] + up + right - up_and_right;
		}
	}

	// The uccs of the particles themselves
	std::vector<int> particle_uccs(num_particles);
	for(int i=0; i<num_particles; i++)
		particle_uccs[i] = ucc[num_particles - r2[i] - 1][r1[i]];

	// Sort the particle uccs from highest to lowest
	std::vector<int> particle_uccs_sorted = particle_uccs;
	std::sort(particle_uccs_sorted.begin(), particle_uccs_sorted.end());
	std::reverse(particle_uccs_sorted.begin(), particle_uccs_sorted.end());

	// Make a ucc threshold (particles on the threshold die too)
	int threshold, threshold_id;
	int count_thresholds_tried = 0;

	threshold_id = num_particles/2;
	threshold_selection:
	threshold = particle_uccs_sorted[threshold_id];
	++count_thresholds_tried;

	// -1 <===> interior
	//  0 <===> boundary
	// +1 <===> exterior
	status.assign(num_particles, 0);
	int num_interior = 0;
	int num_exterior = 0;
	int num_boundary = 0;
	for(int i=0; i<num_particles; i++)
	{
		if(particle_uccs[i] > threshold)
		{
			status[i] = -1;
			++num_interior;
		}
		if(particle_uccs[i] < threshold)
		{
			status[i] = +1;
			++num_exterior;
		}
	}
	num_boundary = num_particles - (num_interior + num_exterior);

	// Handle the case where there are no interior or exterior particles
	if(num_interior == 0 || num_exterior == 0)
	{
		if(count_thresholds_tried == num_particles)
		{
			// Print messages and quit
			std::cout<<"# Iteration "<<(iteration+1)<<"."<<std::endl;
			std::cout<<"# (num_interior, num_boundary, num_exterior) = (";
			std::cout<<num_interior<<", "<<num_boundary<<", "<<num_exterior<<")."<<std::endl;
			std::cout<<"# CANNOT CONTINUE."<<std::endl;
			exit(0);
		}

		if(threshold_id < num_particles/2 && threshold_id != (num_particles-1))
			--threshold_id;
		else if(threshold_id == (num_particles-1))
			threshold_id = num_particles/2 - 1;
		else
			++threshold_id;

		goto threshold_selection;
	}
	else // Standard TwinPeaks
	{
		// Place forbidding rectangles anywhere ucc >= threshold
		for(int i=0; i<num_particles; i++)
		{
			int j = num_particles-1;
			while(j > 0 && ucc[i][j] < threshold)
				j--;

			if(ucc[i][j] >= threshold)
			{
				std::vector<ScalarType> latest{s1[j], s2[num_particles-i-1]};
				prune_rectangles(latest);
				rects.push_front(latest);
			}
		}
	}

	// Select some particles to save in their entirety
	std::vector<bool> save(num_particles, false);
	if(num_interior != 0)
	{
		for(int i=0; i<saves_per_iteration; i++)
		{
			int ii;
			do
			{
				ii = rngs[0].rand_int(num_particles);
			}while(status[ii] != -1);
			save[ii] = true;
		}
	}

	// Save interior particles
	for(int i=0; i<num_particles; i++)
	{
		if(status[i] == -1)
		{
			// Assign prior weight
			double logw = log_prior_mass - log(num_particles - num_boundary);

			// Write it out to an output file
			std::fstream fout;
			if(save[i])
			{
				fout.open("sample.txt", std::ios::out|std::ios::app);
				fout<<logw<<' ';
				for(ScalarType s: scalars[i])
					fout<<s.get_value()<<' ';
				particles[i].write_text(fout);
				fout<<std::endl;
				fout.close();
			}
			fout.open("sample_info.txt", std::ios::out|std::ios::app);
			fout<<logw<<' ';
			for(ScalarType s: scalars[i])
				fout<<s.get_value()<<' ';
			fout<<std::endl;
			fout.close();
		}
	}

	// Reduce remaining prior mass
	log_prior_mass = logdiffexp(log_prior_mass, log_prior_mass + log(num_interior) - log(num_particles - num_boundary));

	// Print messages
	std::cout<<"# Iteration "<<(iteration+1)<<"."<<std::endl;
	std::cout<<"# (num_interior, num_boundary, num_exterior) = (";
	std::cout<<num_interior<<", "<<num_boundary<<", "<<num_exterior<<")."<<std::endl;
	std::cout<<"# Killing "<<(num_interior + num_boundary)<<" particles. "<<std::endl;
	std::cout<<"# "<<rects.size()<<" rectangles. Log(remaining prior mass) = ";
	std::cout<<log_prior_mass<<"."<<std::endl;

	// Replace dead particles
	int num_threads = rngs.size();

	// Backup
	backup_particles = particles;
	backup_scalars = scalars;

	// Assign particles to threads
	std::vector< std::vector<int> > which_particles(num_threads);
	int k=0;
	for(int i=0; i<num_particles; ++i)
	{
		if(status[i] != 1)
			which_particles[(k++)%num_threads].push_back(i);
	}
	// Store acceptance counts
	std::vector<int> accepts(num_threads);

	// Do MCMC to generate new particles
	std::vector<std::thread> threads;
	for(int i=0; i<num_threads; i++)
	{
		threads.emplace_back(std::thread(std::bind(&Sampler<MyModel>::refresh_particles, this, which_particles[i], i, std::ref(accepts[i]))));
	}
	for(int i=0; i<num_threads; i++)
		threads[i].join();
	// Sum acceptance counts
	int accepted = 0;
	for(const int& c: accepts)
		accepted += c;

	std::cout<<"# Accepted "<<accepted<<"/"<<(num_interior + num_boundary)*mcmc_steps<<" (";
	std::cout<<std::fixed<<std::setprecision(1);
	std::cout<<(100.*accepted/((num_interior + num_boundary)*mcmc_steps));
	std::cout<<"%)."<<std::endl<<std::endl;
	std::cout<<std::defaultfloat<<std::setprecision(6);

	iteration++;
	return log_prior_mass;
}


// Refresh all the particles in indices
template<class MyModel>
void Sampler<MyModel>::refresh_particles(const std::vector<int>& indices,
														int which_rng,
														int& accepts)
{
	accepts = 0;
	for(int i: indices)
		accepts += refresh_particle(i, which_rng);
}

template<class MyModel>
int Sampler<MyModel>::refresh_particle(int which, int which_rng)
{
	// Choose a particle to clone
	int copy;
	do
	{
		copy = rngs[which_rng].rand_int(num_particles);
	}
	while(status[copy] != 1);

	// Clone it
	particles[which] = backup_particles[copy];
	scalars[which] = backup_scalars[copy];

	// Do the MCMC
	MyModel proposal;
	std::vector<ScalarType> s_proposal;
	double logH;
	int accepted = 0;
	for(int i=0; i<mcmc_steps; i++)
	{
		proposal = particles[which];
		s_proposal = scalars[which];

		logH = proposal.perturb(rngs[which_rng]);
		for(size_t j=0; j<s_proposal.size(); j++)
			s_proposal[j].set_value(proposal.get_scalars()[j]);
		logH += s_proposal[rngs[which_rng].rand_int(s_proposal.size())].perturb(rngs[which_rng]);

		if(rngs[which_rng].rand() <= exp(logH) && is_okay(s_proposal))
		{
			particles[which] = proposal;
			scalars[which] = s_proposal;
			accepted++;
		}
	}
	return accepted;
}
*/

} // namespace TwinPeaks

