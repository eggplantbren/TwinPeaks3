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

	public:
		MyModel();

		void from_prior();
		double perturb();

		friend std::ostream& operator << (std::ostream& out,
							const MyModel& m);
};

std::ostream& operator << (std::ostream& out, const MyModel& m);

