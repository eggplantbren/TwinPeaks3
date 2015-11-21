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
	assert(num_particles%2 == 0);

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
	// Calculate all upper corner counts
	std::vector<double> uccs(num_particles, 0.5);
	for(int i=0; i<num_particles; i++)
	{
		for(int j=0; j<num_particles; j++)
			uccs[i] += is_in_lower_rectangle(scalars[i], scalars[j]);
	}

	// Argsort the corner counts
	std::vector<size_t> indices = argsort(uccs);
	std::reverse(indices.begin(), indices.end());

	// Mark first half the particles as dying
	std::vector<bool> dying(num_particles, false);
	for(int i=0; i<num_particles/2; i++)
		dying[indices[i]] = true;
	// Continue until ucc changes
	for(int i=num_particles/2; i<num_particles; i++)
	{
		if(uccs[indices[i]] != uccs[indices[num_particles/2-1]])
			break;
		dying[indices[i]] = true;
	}

	// Count number of dying particles
	int num_dying = 0;
	for(const bool d: dying)
		num_dying += static_cast<int>(d);

	// For each dying particle
	for(int i=0; i<num_particles; i++)
	{
		if(dying[i])
		{
			// Append its scalars to the forbidden rectangles
			rects.push_front(scalars[i]);
			prune_rectangles();
		}
	}

	// Mass saved
	double logw_tot = -1E300;

	// Save dying particles
	for(int i=0; i<num_particles; i++)
	{
		// Assign prior weight
		double logw = log_prior_mass - log(num_particles - num_dying);
		logw_tot = logsumexp(logw_tot, logw);

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

//	for(int i=0; i<num_particles; i++)
//	{
//		std::cout<<dying[i]<<' '<<interior[i]<<' '<<uccs[i]<<' ';
//		for(size_t j=0; j<scalars[i].size(); j++)
//			std::cout<<scalars[i][j].get_value()<<' ';
//		std::cout<<std::endl;
//	}
//	exit(0);

	// Reduce remaining prior mass
	log_prior_mass = logdiffexp(log_prior_mass, logw_tot);

	std::cout<<"# Iteration "<<(iteration+1)<<". ";
	std::cout<<"Killing "<<num_dying<<" particles. ";
	std::cout<<"# "<<rects.size()<<" rectangles. Log(remaining prior mass) = ";
	std::cout<<log_prior_mass<<"."<<std::endl;

	int accepted = 0;

	// For each dying particle
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

