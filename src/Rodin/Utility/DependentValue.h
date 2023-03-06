#ifndef RODIN_UTILITY_DEPENDENTVALUE_H
#define RODIN_UTILITY_DEPENDENTVALUE_H

namespace Rodin::Utility
{
  template <class T, T Value_, class ... Args>
  struct DependentValue { static constexpr bool Value = Value_; };
}

#endif


