#include <cassert>
#include "Context.h"

using namespace std;
using namespace TwinPeaks;

Context::Context(double opacity)
:opacity(opacity)
,points{}
{

}

void Context::add_point(const vector<ScalarType>& point)
{
    assert(point.size() == 2);
    points.push_front(point);
}

size_t Context::calculate_ucc(const vector<ScalarType>& point)
{
    assert(point.size() == 2);
    size_t ucc = 0;

    for(auto it=points.begin(); it != points.end(); ++it)
    {
        if(ScalarType::compare(point, *it) == -1)
            ++ucc;
    }
}

