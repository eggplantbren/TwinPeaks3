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

	// Write it out to an output file
	std::fstream fout("sample_info.txt", std::ios::out|std::ios::app);
	const std::vector<double>& scalars = particles[which].get_scalars();
	for(double s: scalars)
		fout<<s<<' ';
	fout<<std::endl;
	fout.close();
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

