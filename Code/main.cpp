#include <iostream>
#include <ctime>
#include <RandomNumberGenerator.h>
#include "Sampler.h"
#include "MyModel.h"

using namespace std;
using namespace DNest3;

int main()
{
	// Initialise RNG and seed with time
	RandomNumberGenerator::initialise_instance();
	RandomNumberGenerator::get_instance().set_seed(time(0));

	Sampler<MyModel> s(50, 1000, 100);
	s.initialise();

	for(int i=0; i<10000; i++)
		s.update();

	return 0;
}

