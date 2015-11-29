#ifndef _Gravity_
#define _Gravity_

#include "RNG.h"
#include <vector>
#include <ostream>

class Gravity
{
	private:
		std::vector<double> x, y, z, vx, vy, vz;

		long double PE, KE, L;
		int staleness;

		std::vector<double> scalars;

		void refresh();
		void increment(int i, int sign);
		void compute_scalars();
	public:
		Gravity();
		void from_prior(RNG& rng);
		double perturb(RNG& rng);

		void write_text(std::ostream& out) const;

		// Getter
		const std::vector<double>& get_scalars() const
		{ return scalars; }
};

#endif

