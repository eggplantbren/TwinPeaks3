#ifndef _Gravity_
#define _Gravity_

#include "Model.h"
#include <vector>
#include <ostream>

class Gravity:public Model
{
	private:
		std::vector<double> x, y, z, vx, vy, vz;

		long double PE, KE, L;
		int staleness;

		void refresh();
		void increment(int i, int sign);
		void compute_scalars();
	public:
		Gravity();
		void from_prior();
		double perturb();

	friend std::ostream& operator << (std::ostream& out, const Gravity& e);
};

#endif

