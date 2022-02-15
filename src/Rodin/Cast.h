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
#include "Rodin/Utility/Overloaded.h"

namespace Rodin
{
   /**
    * @brief Class to help in the casting of objects from a type to another.
    * @tparam From Type to be converted/cast
    * @tparam Mutable Specifies whether the Cast instance will be allowed to
    * mutate the casted object.
    */
   template <class From, bool Mutable = false>
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
   using ImmutableCast = Cast<From, false>;

   /**
    * @internal
    * @brief Deduction guide to aid in the construction of a mutable cast
    * instances, from l-values.
    */
   template <class From>
   Cast(From&) -> Cast<From, false>;

   /**
    * @internal
    * @brief Deduction guide to aid in the construction of a mutable cast
    * instances, from const l-values.
    */
   template <class From>
   Cast(const From&) -> Cast<From, false>;

   /**
    * @brief Specialization representing an immutable cast.
    * @tparam From Type to be converted/cast
    *
    * An immutable cast is that which is not allowed to mutate the underlying
    * casted object. This class should be initialized with l-values of type
    * From.
    */
   template <class From>
   class Cast<From, false>
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
          * @brief Perform cast to type To.
          * @tparam To Target type, to which the object will be casted.
          * @returns Value of type To.
          */
         template <class To>
         To to() const
         {
            return static_cast<To>(from());
         }

      private:
         const From& m_from;
   };

   /**
    * @brief Specialization representing a mutable cast.
    * @tparam From Type to be converted/cast
    *
    * A mutable cast is that which is allowed to mutate the underlying casted
    * object. This class should be initialized with r-values of type From.
    */
   template <class From>
   using MutableCast = Cast<From, true>;

   /**
    * @internal
    * @brief Deduction guide to aid in the construction of a mutable cast
    * instances, from r-values or moved objects.
    */
   template <class From>
   Cast(From&&) -> Cast<From, true>;

   /**
    * @brief Specialization representing an immutable cast.
    * @tparam From Type to be converted/cast
    *
    * An immutable cast is that which is not allowed to mutate the underlying
    * casted object. This class should be initialized with l-values of type
    * From.
    */
   template <class From>
   class Cast<From, true>
   {
      public:
         Cast(From&& from)
            : m_from(std::move(from))
         {}

         const From& from() const
         {
            return m_from;
         }

         From& from()
         {
            return m_from;
         }

         template <class To>
         To to() const
         {
            return static_cast<To>(from());
         }

      private:
         From m_from;
   };
}

#endif
