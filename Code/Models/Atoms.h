#ifndef _Atoms_
#define _Atoms_

#include "Model.h"
#include <vector>

class Atoms:public Model
{
	private:
		// Positions
		std::vector<double> x, y;

		// Energy
		long double PE;

		// From scratch, total
		void calculate_PE();

		// Just a single pair
		double calculate_PE(int i, int j);
		void compute_scalars();

	public:
		Atoms();
		~Atoms();

		// Generate the point from the prior
		void from_prior();

		// Metropolis-Hastings proposals
		double perturb();

		friend std::ostream& operator << (std::ostream& out, const
								Atoms& a);
};

std::ostream& operator << (std::ostream& out, const Atoms& a);

#endif

