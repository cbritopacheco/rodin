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

  template <typename ...>
  struct IsOneOf
  {
    static constexpr bool Value = false;
  };

  template <typename F, typename S, typename ... T>
  struct IsOneOf<F, S, T...> {
    static constexpr bool Value =
      std::is_same<F, S>::Value || IsOneOf<F, T...>::Value;
  };

  template <class T>
  struct IsPlainObject;

  template <class EigenDerived>
  struct IsPlainObject
  {
    static constexpr const bool Value =
      std::is_base_of_v<Eigen::PlainObjectBase<std::decay_t<EigenDerived>>, std::decay_t<EigenDerived>>;
  };

  template <class EigenDerived>
  struct IsPlainObject<Variational::TensorBasis<EigenDerived>>
  {
    static constexpr const bool Value = IsPlainObject<EigenDerived>::Value;
  };
}

#endif
