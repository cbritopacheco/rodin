#ifndef RODIN_UTILITY_WRAP_H
#define RODIN_UTILITY_WRAP_H

#include <utility>

#include "Rodin/Tuple.h"

namespace Rodin::Utility
{
  template <class ...>
  class Wrap;

  template <class ... Ts>
  class Wrap<Tuple<Ts...>>
  {
    public:
      template <template <class> class External>
      using Type = Tuple<External<Ts>...>;
  };

}

#endif

