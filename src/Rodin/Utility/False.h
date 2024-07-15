#ifndef RODIN_UTILITY_FALSE_H
#define RODIN_UTILITY_FALSE_H

namespace Rodin::Utility
{
  template <class...>
  inline constexpr bool False = false;
}

#endif

