#include <iostream>
#include <ctime>
#include <cstdlib>
#include <RandomNumberGenerator.h>
#include "Sampler.h"
#include "Models/SimpleExample.h"

using namespace std;
using namespace DNest3;

int main()
{
	// Initialise RNG and seed with time
	RandomNumberGenerator::initialise_instance();
	RandomNumberGenerator::get_instance().set_seed(time(0));

	Sampler<SimpleExample> sampler(10);
	sampler.initialise();

	for(int i=0; i<1000; i++)
	{
		sampler.explore();
		sampler.refresh();
	}

	return 0;
}

