/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CAST_H
#define RODIN_CAST_H

#include <variant>
#include <utility>

#include "Rodin/Alert.h"
#include "Rodin/Traits.h"
#include "Rodin/Utility/Overloaded.h"

namespace Rodin
{
  template <class From, class DstType>
  class ADLCaster
  {
    public:
      DstType cast(const From& src);
  };

  /**
   * @brief Class to help in the casting of objects from a type to another.
   * @tparam From Type to be converted/cast
   * @tparam Mutable Specifies whether the Cast instance will be allowed to
   * mutate the casted object.
   */
  template <class From>
  class Cast;

  /**
   * @brief Specialization representing an immutable cast.
   * @tparam From Type to be converted/cast
   *
   * An immutable cast is that which is not allowed to mutate the underlying
   * casted object. This class should be initialized with l-values of type
   * From.
   */
  template <class From>
  class Cast
  {
    public:
      /**
       * @brief Constructs a Cast instance with the object to be casted.
       * @param[in] from Object to be casted.
       */
      Cast(const From& from)
        : m_from(from)
      {}

      /**
       * @brief Gets the constant reference to the object which will be
       * casted.
       * @returns Constant reference to the object which will be casted.
       */
      const From& from() const
      {
        return m_from;
      }

      /**
       * @brief Perform cast to type DstType.
       *
       * @tparam DstType Target type to which the object will be casted.
       * @param[in] args List of arguments to be forwarded to the ADLCaster.
       * @returns Value of type To.
       */
      template <class DstType, class ... Args>
      DstType to(Args&&... args) const
      {
        static_assert(std::is_constructible_v<ADLCaster<From, DstType>, Args...>);
        return ADLCaster<From, DstType>(std::forward<Args>(args)...).cast(from());
      }

    private:
      const From& m_from;
  };
}

#endif
