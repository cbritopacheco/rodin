#include "RangeShape.h"

namespace Rodin::Variational
{
  std::ostream& operator<<(std::ostream& os, const RangeShape& obj)
  {
    os << "{" << obj.height() << ", " << obj.width() << "}";
    return os;
  }
}
