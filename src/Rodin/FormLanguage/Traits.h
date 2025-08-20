#ifndef RODIN_FORMLANGUAGE_TRAITS_H
#define RODIN_FORMLANGUAGE_TRAITS_H

#include <type_traits>
#include <boost/type_index.hpp>

#include <Eigen/Core>

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

  template <class T, class Enable = void>
  struct Traits;
}

#endif
