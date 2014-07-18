#ifndef _Sampler_
#define _Sampler_

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

template<class Type>
class Sampler
{
	private:
		std::vector<Type> particles;
		std::vector<double> threshold;
		std::vector<double> direction;

		int mcmc_steps, thin;
		int iteration;

		int find_worst(int which_scalar) const;

	public:
		Sampler(int num_particles, int mcmc_steps, int thin);

		void initialise();
		void update();

};

#include "SamplerImpl.h"

#endif

