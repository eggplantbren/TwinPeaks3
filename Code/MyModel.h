#ifndef _MyModel_
#define _MyModel_

/*
* An object of this class represents a
* point in the parameter space.
*/

#include <vector>
#include <ostream>

class MyModel
{
	private:
		std::vector<double> params;

		void compute_scalars();	
		std::vector<double> scalars;

		std::vector<double> threshold;

	public:
		MyModel();

		void from_prior();
		double perturb();

		// Are all scalars greater than or equal to the threshold
		// values?
		bool is_above(const std::vector<double>& threshold) const;

		// Getter
		const std::vector<double>& get_scalars() const
		{ return scalars; }

		friend std::ostream& operator << (std::ostream& out,
							const MyModel& m);
};

std::ostream& operator << (std::ostream& out, const MyModel& m);

#endif

