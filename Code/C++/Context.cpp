#include "Context.h"

using namespace std;
using namespace TwinPeaks;

Context::Context()
:points{}
{

}

void Context::add_point(const vector<ScalarType>& point)
{
    points.push_front(point);
}

