/**
 * @file IsPlaneObject.h
 * @brief Type trait for identifying Eigen plain object types.
 *
 * This file provides the IsPlainObject type trait, which determines whether
 * a type is an Eigen plain object (Matrix, Array, Vector, etc.) that owns
 * its data, as opposed to expression templates or map types.
 */
#include <type_traits>

#include <Eigen/Core>

#include "Rodin/Types.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Type trait to identify Eigen plain object types.
   * @tparam T Type to check
   * @ingroup RodinFormLanguage
   *
   * This trait determines whether a type is an Eigen plain object type that
   * owns its data storage, such as Matrix, Array, or Vector. Plain objects
   * are distinguished from expression templates and mapped types.
   *
   * ## Type Classification
   * - **Plain objects**: Matrix, Array, Vector (own their data)
   * - **Non-plain**: MatrixBase expressions, Map types, Block expressions
   *
   * ## Usage
   * @code{.cpp}
   * // Check if type is plain object
   * static_assert(IsPlainObject<Eigen::MatrixXd>::Value);
   * static_assert(IsPlainObject<Eigen::VectorXd>::Value);
   * static_assert(!IsPlainObject<Eigen::MatrixXd::RowXpr>::Value);
   * @endcode
   *
   * @note This trait is used internally for object lifetime management in
   * FormLanguage::Base to determine whether objects should be stored or
   * referenced.
   */
  template <class T>
  struct IsPlainObject;

  /**
   * @brief General case: checks if type derives from Eigen::PlainObjectBase.
   * @tparam T Type to check
   * @ingroup RodinFormLanguage
   */
  template <class T>
  struct IsPlainObject
  {
    /**
     * @brief True if T is derived from Eigen::PlainObjectBase<T>, false otherwise.
     */
    static constexpr const bool Value = std::is_base_of_v<Eigen::PlainObjectBase<T>, T>;
  };

  /**
   * @brief Specialization for Eigen::PlainObjectBase types.
   * @tparam Derived Derived type parameter
   * @ingroup RodinFormLanguage
   *
   * This specialization handles the case when T is exactly Eigen::PlainObjectBase<Derived>,
   * ensuring correct trait evaluation for base class references.
   */
  template <class Derived>
  struct IsPlainObject<Eigen::PlainObjectBase<Derived>>
  {
    /**
     * @brief True if Derived is actually derived from PlainObjectBase<Derived>.
     */
    static constexpr const bool Value = std::is_base_of_v<Eigen::PlainObjectBase<Derived>, Derived>;
  };
}
