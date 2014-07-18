#ifndef _SimpleExample_
#define _SimpleExample_

/*
* An object of this class represents a
* point in the parameter space.
*/

#include <vector>
#include <ostream>

class SimpleExample
{
	private:
		std::vector<double> params;

		void compute_scalars();	
		std::vector<double> scalars, tiebreakers;

	public:
		SimpleExample();

		void from_prior();
		double perturb();

		void from_prior_tiebreakers();
		double perturb_tiebreakers();

		// Are all scalars greater than or equal to the threshold
		// values?
		bool is_above(const std::vector<double>& threshold) const;

		// Getter
		const std::vector<double>& get_scalars() const
		{ return scalars; }

		friend std::ostream& operator << (std::ostream& out,
							const SimpleExample& m);
};

std::ostream& operator << (std::ostream& out, const SimpleExample& m);

#endif

