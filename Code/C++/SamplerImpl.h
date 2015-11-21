#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <algorithm>
#include <cassert>
#include "Utils.h"

template<class MyModel>
Sampler<MyModel>::Sampler(const RNG& rng, int num_particles, int mcmc_steps,
							int save_interval)
:rng(rng)
,num_particles(num_particles)
,particles(num_particles)
,scalars(num_particles)
,mcmc_steps(mcmc_steps)
,save_interval(save_interval)
,initialised(false)
,iteration(0)
,log_prior_mass(0.)
{
	assert(num_particles%2 != 0);

	// Open and close outputs file to clear them
	std::fstream fout;
	fout.open("sample.txt", std::ios::out);
	fout.close();
	fout.open("sample_info.txt", std::ios::out);
	fout.close();
}

template<class MyModel>
void Sampler<MyModel>::set_rng_seed(unsigned int seed)
{
	rng.set_seed(seed);
}

template<class MyModel>
void Sampler<MyModel>::initialise()
{
	for(int i=0; i<num_particles; i++)
	{
		particles[i].from_prior(rng);
		const std::vector<double>& s = particles[i].get_scalars();
		scalars[i].clear();
		for(size_t j=0; j<s.size(); j++)
			scalars[i].push_back(ScalarType(s[j]));
		for(size_t j=0; j<scalars[i].size(); j++)
			scalars[i][j].from_prior(rng);
	}
}

template<class MyModel>
void Sampler<MyModel>::prune_rectangles()
{
	// Iterator pointing at first rectangle (most recent)
	auto it0 = rects.begin();

	// Iterator pointing at second rectangle
	auto it=rects.begin();
	it++;

	// Remove redundant rectangles
	for(; it != rects.end(); it++)
	{
		if(is_in_lower_rectangle(*it, *it0))
			it = rects.erase(it);
	}
}

template<class MyModel>
void Sampler<MyModel>::do_iteration()
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
	std::vector< std::vector<int> > empirical_measure(num_particles,
										std::vector<int>(num_particles, 0));
	for(int i=0; i<num_particles; i++)
		empirical_measure[num_particles - r2[i] - 1][r1[i]] += 1;

	// Integrate empirical measure to get (inclusive) upper corner count
	std::vector< std::vector<int> > ucc(num_particles,
										std::vector<int>(num_particles, 0));
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
	int threshold = particle_uccs_sorted[num_particles/2];

	// Assign particles to die
	std::vector<bool> dying(num_particles, false);
	int num_dying = 0;
	for(int i=0; i<num_particles; i++)
	{
		if(particle_uccs[i] >= threshold)
		{
			dying[i] = true;
			num_dying++;
		}
	}

	// Save dying particles
	// Mass saved
	for(int i=0; i<num_particles; i++)
	{
		if(dying[i])
		{
			// Assign prior weight
			double logw = log_prior_mass - log(num_particles);

			// Write it out to an output file
			std::fstream fout;
			if((iteration+1)%save_interval == 0)
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

	// Place forbidding rectangles anywhere ucc >= threshold
	for(int i=0; i<num_particles; i++)
	{
		int j;
		for(j=(num_particles-1); j>=0; j--)
			if(ucc[i][j] >= threshold)
				break;
		if(ucc[i][j] >= threshold)
		{
			rects.push_front({s1[i], s2[j]});
			prune_rectangles();
		}
	}

	// Reduce remaining prior mass
	log_prior_mass = logdiffexp(log_prior_mass, log_prior_mass + log(num_dying) - log(num_particles));

	// Print messages
	std::cout<<"# Iteration "<<(iteration+1)<<". ";
	std::cout<<"Killing "<<num_dying<<" particles. "<<std::endl;
	std::cout<<"# "<<rects.size()<<" rectangles. Log(remaining prior mass) = ";
	std::cout<<log_prior_mass<<"."<<std::endl;

	// Replace dead particles
	int accepted = 0;
	for(int i=0; i<num_particles; ++i)
	{
		if(dying[i])
		{
			// Do MCMC to generate a new particle
			accepted += refresh_particle(i);
		}
	}
	std::cout<<"# Accepted "<<accepted<<"/"<<num_dying*mcmc_steps<<". "<<std::endl<<std::endl;
	iteration++;
}

template<class MyModel>
int Sampler<MyModel>::refresh_particle(int which)
{
	// Choose a particle to clone
	int copy;
	do
	{
		copy = rng.rand_int(num_particles);
	}
	while(!is_okay(scalars[copy]));

	// Clone it
	particles[which] = particles[copy];
	scalars[which] = scalars[copy];

	// Do the MCMC
	MyModel proposal;
	std::vector<ScalarType> s_proposal;
	double logH;
	int accepted = 0;
	for(int i=0; i<mcmc_steps; i++)
	{
		proposal = particles[which];
		s_proposal = scalars[which];

		logH = proposal.perturb(rng);
		for(size_t j=0; j<s_proposal.size(); j++)
			s_proposal[j].set_value(proposal.get_scalars()[j]);
		logH += s_proposal[rng.rand_int(s_proposal.size())].perturb(rng);

		if(rng.rand() <= exp(logH) && is_okay(s_proposal))
		{
			particles[which] = proposal;
			scalars[which] = s_proposal;
			accepted++;
		}
	}
	return accepted;
}

template<class MyModel>
bool Sampler<MyModel>::is_okay(const std::vector<ScalarType>& s)
{
	for(auto it=rects.begin();
				it != rects.end(); it++)
	{
		if(is_in_lower_rectangle(s, *it))
			return false;
	}
	return true;
}

template<class MyModel>
bool Sampler<MyModel>::is_in_lower_rectangle(const std::vector<ScalarType>& s,
											const std::vector<ScalarType>& rect)
{
	for(size_t i=0; i<s.size(); i++)
	{
		if(!(s[i] < rect[i]))
			return false;
	}
	return true;
}

