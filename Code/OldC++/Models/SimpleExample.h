#ifndef _SimpleExample_
#define _SimpleExample_

/*
* An object of this class represents a
* point in the parameter space.
*/

#include <vector>
#include <ostream>
#include "Model.h"

class SimpleExample:public Model
{
	private:
		std::vector<double> params;

		void compute_scalars();	

	public:
		SimpleExample();
		~SimpleExample();

		void from_prior();
		double perturb();

		friend std::ostream& operator << (std::ostream& out,
							const SimpleExample& m);
};

std::ostream& operator << (std::ostream& out, const SimpleExample& m);

#endif

