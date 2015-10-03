#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include <fstream>
#include "Utils.h"

template<class MyModel>
Sampler<MyModel>::Sampler(const RNG& rng, int num_particles, int mcmc_steps)
:rng(rng)
,num_particles(num_particles)
,particles(num_particles)
,mcmc_steps(mcmc_steps)
,initialised(false)
,iteration(0)
{
	// Open and close output file to clear it
	std::fstream fout("sample_info.txt", std::ios::out);
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
	for(MyModel& p: particles)
		p.from_prior(rng);
}

template<class MyModel>
void Sampler<MyModel>::prune_rectangles()
{
	// Iterator pointing at first rectangle (most recent)
	std::list< std::vector<double> >::iterator it0 = rects.begin();

	// Iterator pointing at second rectangle
	std::list< std::vector<double> >::iterator it=rects.begin();
	it++;

	// Remove redundant rectangles
	for(; it != rects.end(); it++)
	{
		if(is_in_rectangle(*it, *it0))
			it = rects.erase(it);
	}
}

template<class MyModel>
void Sampler<MyModel>::do_iteration()
{
	// Calculate how many particles are in the rectangle of each particle
	std::vector<int> counts(num_particles, 0);
	std::vector<int> indices; // Save counts==0 indices
	for(int i=0; i<num_particles; i++)
	{
		for(int j=0; j<num_particles; j++)
		{
			counts[i] += is_in_rectangle(particles[j].get_scalars(),
								particles[i].get_scalars());
		}
		if(counts[i] == 0)
			indices.push_back(i);
	}

	// Choose an index, to become the discarded particle
	int which = indices[rng.rand_int(indices.size())];

	// Append its scalars to the forbidden rectangles
	rects.push_front(particles[which].get_scalars());
	prune_rectangles();

	// Assign prior weight
	double logw = log(1./num_particles) + iteration*log(1. - 1./num_particles);

	// Write it out to an output file
	std::fstream fout("sample_info.txt", std::ios::out|std::ios::app);
	fout<<logw<<' ';
	for(double s: particles[which].get_scalars())
		fout<<s<<' ';
	fout<<std::endl;
	fout.close();

	// Do MCMC to generate a new particle
	refresh_particle(which);
	iteration++;
}

template<class MyModel>
void Sampler<MyModel>::refresh_particle(int which)
{
	// Choose a particle to clone
	int copy;
	do
	{
		copy = rng.rand_int(num_particles);
	}
	while(copy == which && num_particles > 1);

	// Clone it
	particles[which] = particles[copy];

	// Do the MCMC
	MyModel proposal;
	double logH;
	int accepted = 0;
	for(int i=0; i<mcmc_steps; i++)
	{
		proposal = particles[which];
		logH = proposal.perturb(rng);

		if(rng.rand() <= exp(logH) && is_okay(proposal.get_scalars()))
		{
			particles[which] = proposal;
			accepted++;
		}
	}

	std::cout<<"# Iteration "<<(iteration+1)<<". ";
	std::cout<<"Accepted "<<accepted<<"/"<<mcmc_steps<<". ";
	std::cout<<rects.size()<<" rectangles."<<std::endl;
}

template<class MyModel>
bool Sampler<MyModel>::is_okay(const std::vector<double>& s)
{
	for(std::list< std::vector<double> >::iterator it=rects.begin();
				it != rects.end(); it++)
	{
		if(is_in_rectangle(s, *it))
			return false;
	}
	return true;
}

template<class MyModel>
bool Sampler<MyModel>::is_in_rectangle(const std::vector<double>& s,
										const std::vector<double>& rect)
{
	for(size_t i=0; i<s.size(); i++)
	{
		if(s[i] >= rect[i])
			return false;
	}
	return true;
}

