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
        std::list< std::vector<ScalarType> > rectangles;
        std::list< double > opacities;

    public:
        // Points starts empty
        Context();

        // Add a new rectangle
        void add_rectangle(const std::vector<ScalarType>& latest,
                           double opacity);

        // Calculate the log-prob of a point according to this context
        double log_prob(const std::vector<ScalarType>& point) const;

        size_t get_num_rectangles() const
        { return rectangles.size(); }
};

} // namespace TwinPeaks

#endif

