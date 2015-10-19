#include "ScalarType.h"

ScalarType::ScalarType()
{
}

ScalarType::ScalarType(double logL, double tieBreaker)
:logL(logL)
,tieBreaker(tieBreaker)
{
}

bool operator < (const ScalarType& l1, const ScalarType& l2)
{
	bool result = false;
	if(l1.logL < l2.logL)
		result = true;
	else if(l1.logL == l2.logL && l1.tieBreaker < l2.tieBreaker)
		result = true;
	return result;
}

