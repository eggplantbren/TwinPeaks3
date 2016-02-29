#ifndef TwinPeaks_Rectangle
#define TwinPeaks_Rectangle

#include <vector>
#include "ScalarType.h"

namespace TwinPeaks
{

/*
*   An object of this class defines a "forbidden rectangle" in the
*   L1-L2 space. There is also a degree of forbiddenness, represented
*   by the 'opacity' variable \in [0, 1].
*/
class Rectangle
{
    private:
        std::vector<ScalarType> scalars;
        double opacity;

    public:
        Rectangle();
        Rectangle(const std::vector<ScalarType>& scalars, double opacity);

        // Getters
        const std::vector<ScalarType>& get_scalars() const
        { return scalars; }
        double get_opacity() const
        { return opacity; }

};

} // namespace TwinPeaks

#endif

