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
    // Check everything's ok with the input
    assert(latest.size() == 2 && opacity >= 0.0 && opacity <= 1.0);

    /* Remove redundant rectangles */
    if(opacity == 1.0)
    {
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
    }

    // Add
    rectangles.push_front(latest);
    opacities.push_front(opacity);
}

// Add a whole list of opaque rectangles
void Context::add_opaque_rectangles(vector< vector<ScalarType> >& latest)
{
    // Remove redundant rectangles from 'latest' before doing anything else
    std::vector<bool> redundant(latest.size(), false);
    short result;
    for(size_t i=0; i<latest.size(); ++i)
    {
        for(size_t j=(i+1); j<latest.size(); ++j)
        {
            result = ScalarType::compare(latest[i], latest[j]);
            if(result == -1)
                redundant[i] = true;
            else if(result == 1)
                redundant[j] = true;
        }
    }

    auto it1 = latest.begin();
    auto it2 = redundant.begin();

    while(it1 != latest.end())
    {
        if(*it2) // If redundant
        {
            it1 = latest.erase(it1);
            it2 = redundant.erase(it2);
        }
        ++it1;
        ++it2;
    }

    // Add the rectangles
    for(it1 = latest.begin(); it1 != latest.end(); ++it1)
        add_rectangle(*it1, 1.0);
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

