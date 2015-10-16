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

	// Create a sampler
	Sampler<SimpleExample> sampler(rng, 1000, 10000, 1000);
	sampler.initialise();

	// Do NS indefinitely
	while(true)
		sampler.do_iteration();

	return 0;
}

