#ifndef TwinPeaks_ScalarType
#define TwinPeaks_ScalarType

#include <istream>
#include <ostream>
#include <vector>
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

        // Print to stream
        void print(std::ostream& out) const;

        // Compare to another ScalarType object
        // Return 1 if first argument is greater, -1 if less, 0 if incomparable
        static short compare(const ScalarType& s1, const ScalarType& s2);

        // Compare vectors of ScalarTypes
        // Return 1 if first argument is greater, -1 if less, 0 if incomparable
        static short compare(const std::vector<ScalarType>& s1,
                             const std::vector<ScalarType>& s2);
};

} // namespace TwinPeaks

// Less-than operator (uses compare)
bool operator < (const TwinPeaks::ScalarType& s1,
                    const TwinPeaks::ScalarType& s2);

#endif

