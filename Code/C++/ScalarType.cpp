#include "ScalarType.h"
#include "Utils.h"
#include <cassert>

namespace TwinPeaks
{

ScalarType::ScalarType()
{
}

ScalarType::ScalarType(double value)
:value(value)
{
}

ScalarType::ScalarType(double value, double tb)
:value(value)
,tiebreaker(tb)
{
	assert(tiebreaker >= 0. && tiebreaker <= 1.);
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

short ScalarType::compare(const ScalarType& other) const
{
    // Check for <
    if(value < other.value)
        return -1;
    if(value == other.value && tiebreaker < other.tiebreaker)
        return -1;

    // Check for >
    if(value > other.value)
        return +1;
    if(value == other.value && tiebreaker > other.tiebreaker)
        return +1;

    // Otherwise, they aren't comparable
    return 0;
}

} // namespace TwinPeaks

