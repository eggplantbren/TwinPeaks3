#include <cassert>
#include <limits>
#include "Context.h"

using namespace std;
using namespace TwinPeaks;

Context::Context()
:rectangles{}
,opacities{}
{

}

void Context::add_rectangle(const vector<ScalarType>& latest, double opacity)
{
    assert(latest.size() == 2 && opacity >= 0.0 && opacity <= 1.0);

    /* Remove redundant rectangles */
    // Two iterators
    auto it1 = rectangles.begin();
    auto it2 = opacities.begin();

    while(it1 != rectangles.end())
	{
		if(ScalarType::compare(latest, *it1) == 1)
        {
			it1 = rectangles.erase(it1);
            it2 = opacities.erase(it2);
        }
        ++it1;
        ++it2;
    }

    rectangles.push_front(latest);
}

double Context::log_prob(const vector<ScalarType>& point) const
{
    double logp = 0.0;

    // Two iterators
    auto it1 = rectangles.begin();
    auto it2 = opacities.begin();

    while(it1 != rectangles.end())
	{
        if(ScalarType::compare(point, *it1) == -1)
        {
            if(*it2 >= 1.0)
                return -std::numeric_limits<double>::max();
            else
                logp += log(1.0 - *it2);
        }

        ++it1;
        ++it2;
    }

    return logp;
}

