#ifndef RODIN_FORMLANGUAGE_TRAITS_H
#define RODIN_FORMLANGUAGE_TRAITS_H

#include <type_traits>
#include <boost/type_index.hpp>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "Rodin/Types.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::FormLanguage
{
  /**
   * @defgroup TraitsSpecializations Traits Template Specializations
   * @brief Template specializations of the Traits class.
   *
   * @see Traits
   */

  template <class ... Args>
  struct Traits;
}

#endif
