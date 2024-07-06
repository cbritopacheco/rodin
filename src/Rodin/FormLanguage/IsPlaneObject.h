#include <type_traits>

#include <Eigen/Core>

#include <unsupported/Eigen/CXX11/Tensor>

#include "Rodin/Types.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::FormLanguage
{
  template <class T>
  struct IsPlainObject;

  template <class T>
  struct IsPlainObject
  {
    static constexpr const bool Value = std::is_base_of_v<Eigen::PlainObjectBase<T>, T>;
  };

  template <class Derived>
  struct IsPlainObject<Eigen::PlainObjectBase<Derived>>
  {
    static constexpr const bool Value = std::is_base_of_v<Eigen::PlainObjectBase<Derived>, Derived>;
  };
}
