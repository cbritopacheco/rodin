#ifndef RODIN_FORMLANGUAGE_TRAITS_H
#define RODIN_FORMLANGUAGE_TRAITS_H

#include "Rodin/Variational/ForwardDecls.h"

#include <type_traits>
#include <Eigen/Core>

namespace Rodin::FormLanguage
{
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
