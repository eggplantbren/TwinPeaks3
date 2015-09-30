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
	RNG rng;

	Sampler<SimpleExample> sampler(rng, 100, 5000);
	sampler.initialise();


	sampler.do_iteration();

	return 0;
}

