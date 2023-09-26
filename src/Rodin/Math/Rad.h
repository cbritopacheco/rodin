#ifndef RODIN_MATH_RAD_H
#define RODIN_MATH_RAD_H

#include "Rodin/Types.h"
#include "Unit.h"

namespace Rodin::Math
{
  /**
   * @brief Represents an angle in radians.
   */
  class Rad : public Unit<Rad, Scalar>
  {
    public:
      using Parent = Unit<Rad, Scalar>;
      using Parent::Parent;
  };
}

#endif
