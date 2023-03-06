#ifndef RODIN_UTILITY_DEPENDENTFALSE_H
#define RODIN_UTILITY_DEPENDENTFALSE_H

namespace Rodin::Utility
{
  template <class ... Args>
  struct DependentFalse { static constexpr bool Value = false; };
}

#endif


