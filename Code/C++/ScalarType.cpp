#include "ScalarType.h"
#include "Utils.h"

ScalarType::ScalarType()
{
}

ScalarType::ScalarType(double value)
:value(value)
{
}

void ScalarType::from_prior(RNG& rng)
{
	tiebreaker = rng.rand();
}

double ScalarType::perturb(RNG& rng)
{
	tiebreaker += rng.randh();
	wrap(tiebreaker, 0., 1.);
	return 0.;
}

bool operator < (const ScalarType& s1, const ScalarType& s2)
{
	bool result = false;
	if(s1.value < s2.value)
		result = true;
	else if(s1.value == s2.value && s1.tiebreaker < s2.tiebreaker)
		result = true;
	return result;
}

