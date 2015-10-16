#ifndef _Atoms_
#define _Atoms_

#include "Model.h"
#include "RNG.h"
#include <vector>

class Atoms
{
	private:
		// Positions
		std::vector<double> x, y, z;

		std::vector< std::vector<double> > terms;

		std::vector<double> scalars;

		// Energy
		double PE;

		// From scratch, total
		void calculate_PE();

		// Just a single pair
		void calculate_PE(int i, int j);
		void compute_scalars();

	public:
		Atoms();
		~Atoms();

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

