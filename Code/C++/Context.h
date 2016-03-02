#ifndef TwinPeaks_Context
#define TwinPeaks_Context

#include <list>
#include <vector>
#include "ScalarType.h"

namespace TwinPeaks
{

/*
* A context is a collection of points, with respect to which
* an upper corner count can be defined.
*/

class Context
{
    private:
        std::list< std::vector<ScalarType> > points;

    public:
        Context();

        // Add a point to the context
        void add_point(const std::vector<ScalarType>& point);

        // Calculate the ucc of a point
        size_t calculate_ucc(const std::vector<ScalarType>& point) const;
};

} // namespace TwinPeaks

#endif

