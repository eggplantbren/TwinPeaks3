/*
* Creates a single sampler and does a run. No-frills.
*/

#include <iostream>
#include <ctime>
#include "RNG.h"
#include "Sampler.h"
#include "Models/SimpleExample.h"

using namespace TwinPeaks;
using namespace std;

int main(int argc, char** argv)
{
	// Number of threads is 1 or a number given on the command line
	int num_threads = (argc==1)?(1):(atoi(argv[1]));

	// Make (num_threads) RNGs: one for the main process
	// and num_threads for use whenever we send out threads
	vector<RNG> rngs(num_threads, RNG());
	auto seed = time(0);
	for(RNG& r: rngs)
		r.set_seed(++seed);

	constexpr int num_particles = 1000;
	constexpr int num_mcmc_steps = 1000;

	// Create a sampler
	Sampler<SimpleExample> sampler(num_particles, num_mcmc_steps, rngs);
	sampler.initialise();

    for(int i=0; i<1000000; ++i)
        sampler.do_iteration();

	return 0;
}

