#include "Rectangle.h"

namespace TwinPeaks
{

Rectangle::Rectangle()
{
}

Rectangle::Rectangle(const std::vector<ScalarType>& scalars, double opacity)
:scalars(scalars)
,opacity(opacity)
{

}

} // namespace TwinPeaks

