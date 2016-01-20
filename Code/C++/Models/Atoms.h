#ifndef _Atoms_
#define _Atoms_

#include "RNG.h"
#include <ostream>
#include <vector>

class Atoms
{
	private:
		static constexpr int num_atoms = 20;
		static constexpr double L = 10.;

		// Positions
		std::vector<double> x, y, z;
		std::vector< std::vector<double> > terms1;
		std::vector< std::vector<double> > terms2;

		std::vector<double> scalars;

		// Energy
		double PE1, PE2;

		// From scratch, total
		void calculate_PE();

		// Just a single pair
		void calculate_PE(int i, int j);
		void compute_scalars();

	public:
		Atoms();

		// Generate the point from the prior
		void from_prior(RNG& rng);

		// Metropolis-Hastings proposals
		double perturb(RNG& rng);

		void write_text(std::ostream& out) const;

		// Getter
		const std::vector<double>& get_scalars() const
		{ return scalars; }
};

#endif

