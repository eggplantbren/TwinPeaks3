#ifndef _Potts_
#define _Potts_

#include <vector>
#include <ostream>
#include "RNG.h"

class Potts
{
	private:
		static const int num_colors = 5;

		std::vector< std::vector<int> > x;
		std::vector<double> scalars;

		// So we can just update on the fly
		int score;
		int score2;

		void compute_score();
		void compute_scalars();
	public:
		Potts();
		void from_prior(RNG& rng);
		double perturb(RNG& rng);

		void write_text(std::ostream& out) const;
		const std::vector<double>& get_scalars() const { return scalars; }
};

#endif

