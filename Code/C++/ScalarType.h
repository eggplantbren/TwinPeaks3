#ifndef _ScalarType_
#define _ScalarType_

#include <istream>
#include <ostream>

class ScalarType
{
	private:
		double logL;
		double tieBreaker;

	public:
		ScalarType();
		ScalarType(double logL, double tieBreaker);

	friend bool operator < (const ScalarType& l1, const ScalarType& l2);
};

bool operator < (const ScalarType& l1, const ScalarType& l2);

#endif

