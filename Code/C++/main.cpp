/*
* Creates a single sampler and does a run. No-frills.
*/

#include <iostream>
#include <ctime>
#include "RNG.h"
#include "Sampler.h"
#include "Models/SimpleExample.h"

using namespace std;

int main()
{
	// Make an RNG
	RNG rng;
	rng.set_seed(time(0));

	constexpr int num_particles = 3001;
	constexpr int num_mcmc_steps = 1000;
	constexpr double depth = 1000.;
	constexpr int steps = depth/log(2.);

	// Create a sampler
	Sampler<SimpleExample> sampler(rng, num_particles, num_mcmc_steps,
											num_particles);
	sampler.initialise();

	// Do NS indefinitely
	for(int i=0; i<steps; ++i)
		sampler.do_iteration();

	return 0;
}

