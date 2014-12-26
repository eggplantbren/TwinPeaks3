#ifndef _Sampler_
#define _Sampler_

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Utils.h"

template<class Type>
class Sampler
{
	private:
		int num_particles;
		std::vector<Type> particles;

		std::vector< std::vector<double> > thresholds;

		bool is_below(const Type& particle, int threshold) const;
		bool is_below(const std::vector<double>& s1,
				const std::vector<double>& s2) const;

		void create_threshold(const std::vector< std::vector<double> >&
						keep);

	public:
		Sampler(int num_particles);

		void initialise();
		void explore();

		static bool greater(const std::vector<double>& s1,
					const std::vector<double>& s2);
};

#include "SamplerImpl.h"

#endif

