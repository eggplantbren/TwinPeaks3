#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <algorithm>
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
	std::list< std::vector<ScalarType> >::iterator it0 = rects.begin();

	// Iterator pointing at second rectangle
	std::list< std::vector<ScalarType> >::iterator it=rects.begin();
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
	// Calculate counts
	std::vector<int> lccs(num_particles, 0);
	for(int i=0; i<num_particles; i++)
	{
		for(int j=0; j<num_particles; j++)
		{
			if(is_in_lower_rectangle(scalars[j], scalars[i]))
				lccs[i]++;
		}
		if(lccs[i] == 0)
			lccs[i] = num_particles;
	}
	int min_lcc = *min_element(lccs.begin(), lccs.end());
	if(min_lcc == num_particles)
	{
		std::cerr<<"# All particles have corner counts of zero. Aborting."<<std::endl;
		exit(1);
	}

	std::vector<int> usable_indices;
	for(int i=0; i<num_particles; i++)
		if(lccs[i] == min_lcc)
			usable_indices.push_back(i);

	// Choose one at random
	int choice = usable_indices[rng.rand_int(usable_indices.size())];

	// Find min of scalar 2 among particles in the rectangle
	int which_scalar = rng.rand_int(2);
	std::vector<ScalarType> ss;
	std::vector<int> indices;
	for(int i=0; i<num_particles; i++)
	{
		if(is_in_lower_rectangle(scalars[i], scalars[choice]))
		{
			ss.push_back(scalars[i][which_scalar]);
			indices.push_back(i);
		}
	}

	ScalarType ss_min = *min_element(ss.begin(), ss.end());
	int which = 0;
	for(size_t i=0; i<ss.size(); i++)
		if(ss[i] == ss_min)
			which = indices[i];

	std::vector<ScalarType> forbid(2);
	forbid = scalars[choice];
	forbid[which_scalar] = scalars[which][which_scalar];

	// Append its scalars to the forbidden rectangles
	rects.push_front(forbid);
	prune_rectangles();

	// Assign prior weight
	double logw = log((double)1./(num_particles + 1)) + log_prior_mass;
	log_prior_mass = logdiffexp(log_prior_mass, logw);

	// Write it out to an output file
	std::fstream fout;
	if((iteration+1)%save_interval == 0)
	{
		fout.open("sample.txt", std::ios::out|std::ios::app);
		fout<<logw<<' ';
		for(ScalarType s: scalars[which])
			fout<<s.get_value()<<' ';
		particles[which].write_text(fout);
		fout<<std::endl;
		fout.close();
	}
	fout.open("sample_info.txt", std::ios::out|std::ios::app);
	fout<<logw<<' ';
	for(ScalarType s: scalars[which])
		fout<<s.get_value()<<' ';
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
		logH += s_proposal[rng.rand_int(s_proposal.size())].perturb(rng);

		if(rng.rand() <= exp(logH) && is_okay(s_proposal))
		{
			particles[which] = proposal;
			scalars[which] = s_proposal;
			accepted++;
		}
	}

	std::cout<<"# Iteration "<<(iteration+1)<<". ";
	std::cout<<rects.size()<<" rectangles. Log(remaining prior mass) = ";
	std::cout<<log_prior_mass<<"."<<std::endl;
	std::cout<<"Accepted "<<accepted<<"/"<<mcmc_steps<<". "<<std::endl<<std::endl;
}

template<class MyModel>
bool Sampler<MyModel>::is_okay(const std::vector<ScalarType>& s)
{
	for(std::list< std::vector<ScalarType> >::iterator it=rects.begin();
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

