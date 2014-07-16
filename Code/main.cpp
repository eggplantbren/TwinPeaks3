#include <iostream>
#include <RandomNumberGenerator.h>
#include "Sampler.h"
#include "MyModel.h"

using namespace std;
using namespace DNest3;

int main()
{
	RandomNumberGenerator::initialise_instance();

	Sampler<MyModel> s(100, 1000);
	s.initialise();

	for(int i=0; i<10000; i++)
		s.update();

	return 0;
}

