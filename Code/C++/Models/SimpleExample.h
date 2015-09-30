#ifndef _SimpleExample_
#define _SimpleExample_

#include <vector>
#include <iostream>
#include "RNG.h"

class SimpleExample
{
	private:
		std::vector<double> params;
		std::vector<double> scalars;

		void compute_scalars();

	public:
		SimpleExample();

		// Generate the point from the prior
		void from_prior(RNG& rng);

		// Metropolis-Hastings proposals
		double perturb(RNG& rng);

		// Likelihood function
		double log_likelihood() const;

		// Write to stream (text format)
		void write_text(std::ostream& out) const;

		// Getter
		const std::vector<double>& get_scalars() const
		{ return scalars; }
};

#endif

