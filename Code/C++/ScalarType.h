#ifndef _ScalarType_
#define _ScalarType_

#include <istream>
#include <ostream>
#include "RNG.h"

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

		bool operator == (const ScalarType& other) const
		{ return (value == other.value && tiebreaker == other.tiebreaker); }

	friend bool operator < (const ScalarType& l1, const ScalarType& l2);
	friend bool operator <= (const ScalarType& l1, const ScalarType& l2);
};

bool operator < (const ScalarType& l1, const ScalarType& l2);
bool operator <= (const ScalarType& l1, const ScalarType& l2);

#endif

