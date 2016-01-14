#ifndef _CambridgeLJ_
#define _CambridgeLJ_

#include "RNG.h"
#include <vector>
#include <valarray>
#include <ostream>

class CambridgeLJ
{
	private:
		int N;
		std::valarray<double> s;
		std::vector<double> scalars;

		void compute_scalars();

	public:
		CambridgeLJ();

		// Generate the point from the prior
		void from_prior(RNG& rng);

		// Metropolis-Hastings proposals
		double perturb(RNG& rng);

		void write_text(std::ostream& out) const;

		// Getter
		const std::vector<double>& get_scalars() const
		{ return scalars; }

		CambridgeLJ& operator = (const CambridgeLJ& other);
};

#endif

