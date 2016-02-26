#ifndef TwinPeaks_ScalarType
#define TwinPeaks_ScalarType

#include <istream>
#include <ostream>
#include "RNG.h"

namespace TwinPeaks
{

class ScalarType
{
	private:
		double value;
		double tiebreaker;

	public:
		ScalarType();
		ScalarType(double val);
		ScalarType(double val, double tb);

		// Randomise the tiebreaker
		void from_prior(RNG& rng);
		double perturb(RNG& rng);

		// Getter and setter
		double get_value() const { return value; }
		void set_value(double val) { value = val; }

        // Compare to another ScalarType object
        // Return 1 if this is greater, -1 if less, 0 if incomparable
        short compare(const ScalarType& other) const;
};

} // namespace TwinPeaks

#endif

