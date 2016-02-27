#include "ScalarType.h"
#include "Utils.h"
#include <cassert>

using namespace std;

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

short ScalarType::compare(const ScalarType& s1, const ScalarType& s2)
{
    // Check for <
    if(s1.value < s2.value)
        return -1;
    if(s1.value == s2.value && s1.tiebreaker < s2.tiebreaker)
        return -1;

    // Check for >
    if(s1.value > s2.value)
        return +1;
    if(s1.value == s2.value && s1.tiebreaker > s2.tiebreaker)
        return +1;

    // Otherwise, they aren't comparable
    return 0;
}

short ScalarType::compare(const vector<ScalarType>& s1,
                          const vector<ScalarType>& s2)
{
    assert(s1.size() == s2.size());

    // Compare element-wise
    std::vector<short> comparisons(s1.size(), 0);
    for(size_t i=0; i<s1.size(); ++i)
    {
        comparisons[i] = compare(s1[i], s2[i]);
        if(i > 0 && comparisons[i] != comparisons[i-1])
            return 0;
    }
    return comparisons[0];
}

} // namespace TwinPeaks

