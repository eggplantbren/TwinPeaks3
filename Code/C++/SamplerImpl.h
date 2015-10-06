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
,ucc_tiebreakers(num_particles)
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
	for(double& x: ucc_tiebreakers)
		x = rng.rand();
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
		if(is_in_lower_rectangle(*it, *it0))
			it = rects.erase(it);
	}
}

template<class MyModel>
int Sampler<MyModel>::upper_corner_count(const MyModel& particle) const
{
	int count = 0;
	for(int j=0; j<num_particles; j++)
	{
		count += is_in_upper_rectangle(particles[j].get_scalars(),
														particle.get_scalars());
	}
	return count;
}

template<class MyModel>
void Sampler<MyModel>::do_iteration()
{
	// Calculate *upper* corner counts
	std::vector<int> corner_counts(num_particles, 0);
	for(int i=0; i<num_particles; i++)
		corner_counts[i] = upper_corner_count(particles[i]);

	// Find the particle with the highest augmented corner count
	int which = 0;
	for(int i=1; i<num_particles; i++)
		if((corner_counts[i] > corner_counts[which]) ||
			((corner_counts[i] == corner_counts[which]) &&
			(ucc_tiebreakers[i] > ucc_tiebreakers[which])))
			which = i;

	// Count number of ties
	int ties = 0;
	for(int i=0; i<num_particles; i++)
		if(i != which)
			if(corner_counts[i] == corner_counts[which])
				ties++;

	// Append its scalars to the forbidden rectangles
	if(ties == 0)
	{
		rects.push_front(particles[which].get_scalars());
		prune_rectangles();
	}

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
	if(ties == 0)
		refresh_particle(which, -1, 1.);
	else
		refresh_particle(which, corner_counts[which], ucc_tiebreakers[which]);
	iteration++;
}

// If ucc_threshold is -1, ignore it.
template<class MyModel>
void Sampler<MyModel>::refresh_particle(int which, int ucc_threshold,
														double tb_threshold)
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
	ucc_tiebreakers[which] = ucc_tiebreakers[copy];

	// Do the MCMC
	MyModel proposal;
	int ucc_proposal = 0;
	double ucc_tiebreaker_proposal;
	double logH;
	int accepted = 0;
	for(int i=0; i<mcmc_steps; i++)
	{
		// Perturb the particle		
		proposal = particles[which];
		logH = proposal.perturb(rng);

		// Perturb the ucc_tiebreaker
		ucc_tiebreaker_proposal = ucc_tiebreakers[which] + rng.randh();
		wrap(ucc_tiebreaker_proposal, 0., 1.);

		if(ucc_threshold != -1)
			ucc_proposal = upper_corner_count(proposal);

		if(rng.rand() <= exp(logH) && is_okay(proposal.get_scalars()))
		{
			// Extra criteria apply if ucc_threshold != -1
			if(ucc_threshold == -1 || (ucc_proposal < ucc_threshold)
				|| ((ucc_proposal == ucc_threshold) &&
					(ucc_tiebreaker_proposal < tb_threshold)))
			{
				particles[which] = proposal;
				accepted++;
			}
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
		if(is_in_lower_rectangle(s, *it))
			return false;
	}
	return true;
}

template<class MyModel>
bool Sampler<MyModel>::is_in_lower_rectangle(const std::vector<double>& s,
												const std::vector<double>& rect)
{
	for(size_t i=0; i<s.size(); i++)
	{
		if(s[i] >= rect[i])
			return false;
	}
	return true;
}

template<class MyModel>
bool Sampler<MyModel>::is_in_upper_rectangle(const std::vector<double>& s,
												const std::vector<double>& rect)
{
	for(size_t i=0; i<s.size(); i++)
	{
		if(s[i] <= rect[i])
			return false;
	}
	return true;
}

