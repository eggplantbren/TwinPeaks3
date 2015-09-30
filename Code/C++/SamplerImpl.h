#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include "Utils.h"

template<class MyModel>
Sampler<Model>::Sampler(const RNG& rng, int num_particles, int mcmc_steps)
:rng(rng)
,num_particles(num_particles)
,particles(num_particles)
,mcmc_steps(mcmc_steps)
,initialised(false)
,iteration(0)
{

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
		p.from_prior();
}

