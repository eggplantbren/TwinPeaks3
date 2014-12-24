#ifndef _Sampler_
#define _Sampler_

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

template<class Type>
class Sampler
{
	private:
		int num_particles;
		std::vector<Type> particles;

	public:
		Sampler(int num_particles);

		void initialise();
};

#include "SamplerImpl.h"

#endif

