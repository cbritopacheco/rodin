#ifndef RODIN_TYPES_H
#define RODIN_TYPES_H

#include <cstddef>
#include <Eigen/Core>

#include "Rodin/Configure.h"

namespace Rodin
{
  using Boolean = bool;
  using Scalar = double;
  using Index = std::size_t;
  using Attribute = std::size_t;

#if __cpp_size_t_suffix < 202011L
  constexpr
  std::size_t operator "" _UZ (unsigned long long x)
  {
    return x;
  }
#endif
}

#endif
