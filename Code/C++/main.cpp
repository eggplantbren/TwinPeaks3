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

//	ImageEntropy::load_data();
//	Sampler<ImageEntropy> sampler(8, 2, 1000, 0.5, 10);

	Sampler<SimpleExample> sampler(8, 1000, 5000, 0.001, 1000);
	sampler.initialise();

	while(true)
	{
		sampler.explore();
		sampler.refresh();
	}

	return 0;
}

