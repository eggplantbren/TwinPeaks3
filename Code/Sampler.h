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
		const int num_threads, num_particles, steps, thin;
		const double peel_factor;

		int iterations;
		std::vector<Type> particles;

		std::vector< std::vector<double> > thresholds;
		std::vector< std::vector<double> > thresholds_tiebreakers;
		double log_prior_mass;

		// How many of the thresholds is the particle below?
		int badness(const Type& particle) const;
		// Is s1 below s2?
		bool is_below(const std::vector<double>& s1,
				const std::vector<double>& s2,
				const std::vector<double>& tb1,
				const std::vector<double>& tb2) const;

		void create_threshold(const std::vector< std::vector<double> >&
						keep,
					const std::vector< std::vector<double> >&
						tiebreakers);
		void remove_redundant_thresholds();

	public:
		Sampler(int num_threads, int num_particles, int steps, double peel_factor,
				int thin);

		void initialise();
		void explore();
		void refresh();
};

#include "SamplerImpl.h"

#endif

