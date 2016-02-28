#ifndef TwinPeaks_Rectangle
#define TwinPeaks_Rectangle

#include <vector>
#include "ScalarType.h"

namespace TwinPeaks
{

/*
*   An object of this class defines a "forbidden rectangle" in the
*   L1-L2 space. There is also a degree of forbiddenness.
*/
class Rectangle
{
    private:
        std::vector<ScalarType> scalars;
        double transparency;

    public:
        Rectangle(const ScalarType& s1, const ScalarType& s2,
                        double transparency);

};

} // namespace TwinPeaks

#endif

