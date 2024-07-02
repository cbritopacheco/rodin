#ifndef RODIN_MATH_RAD_H
#define RODIN_MATH_RAD_H

#include "Rodin/Types.h"
#include "Unit.h"

namespace Rodin::Math
{
  /**
   * @brief Represents an angle in radians.
   */
  class Rad : public Unit<Rad, Real>
  {
    public:
      using Parent = Unit<Rad, Real>;
      using Parent::Parent;
  };
}

#endif
