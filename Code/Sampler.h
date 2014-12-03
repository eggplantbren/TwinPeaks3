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
		std::vector<Type> particles;
		std::vector< std::vector<double> > threshold;
		std::vector<double> direction;

		int num_steps, mcmc_steps, thin;
		int iteration;

		int find_worst(int which_scalar) const;
		void update();

	public:
		Sampler(int num_particles, int num_steps, int mcmc_steps, int thin);

		void initialise();
		void run();
};

#include "SamplerImpl.h"

#endif

