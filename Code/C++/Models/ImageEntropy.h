#ifndef TwinPeaks3_ImageEntropy
#define TwinPeaks3_ImageEntropy

#include <vector>
#include <iostream>
#include "RNG.h"
#include "PSF.h"

class ImageEntropy
{
	private:
		std::vector<double> params;
		std::vector<double> scalars;

		void compute_scalars();

	public:
		ImageEntropy();

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

