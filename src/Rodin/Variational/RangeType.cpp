/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "RangeType.h"

namespace Rodin::Variational
{
  std::ostream& operator<<(std::ostream& os, const RangeType& obj)
  {
    switch (obj)
    {
      case RangeType::Boolean:
      {
        os << "Boolean";
        break;
      }
      case RangeType::Integer:
      {
        os << "Boolean";
        break;
      }
      case RangeType::Real:
      {
        os << "Real";
        break;
      }
      case RangeType::Vector:
      {
        os << "Vector";
        break;
      }
      case RangeType::Matrix:
      {
        os << "Matrix";
        break;
      }
    }
    return os;
  }
}
