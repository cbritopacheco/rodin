/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_FORMLANGUAGE_BASE_H
#define RODIN_FORMLANGUAGE_BASE_H

#include <atomic>
#include <vector>
#include <cassert>
#include <variant>
#include <typeinfo>

#include "Rodin/Types.h"
#include "Rodin/Copyable.h"
#include "Rodin/Identifiable.h"
#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Traits.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "Traits.h"
#include "IsPlaneObject.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Base class for all objects in Rodin's FormLanguage system.
   *
   * This class serves as the foundational base for all form language objects,
   * providing essential services such as unique identification, object lifetime
   * management, and polymorphic copying capabilities. All form language expressions,
   * functions, and operators derive from this base class.
   *
   * @warning The FormLanguage::Base class is not thread safe. Only one thread
   * should access the methods of a FormLanguage object at a time.
   *
   * ## Key Features
   * - **Unique Identification**: Each instance receives a unique UUID for tracking
 * - **Object Management**: Value materialization for expression objects
   * - **Polymorphic Operations**: Support for copying and cloning operations
   * - **Type Safety**: Template-based object storage with type validation
   */
  class Base : public Copyable, public Identifiable
  {
    public:
      /**
       * @brief Constructor.
       */
      Base() = default;

      /**
       * @brief Copy constructor.
       */
      Base(const Base& other)
        : Copyable(other),
          Identifiable(other)
      {}

      /**
       * @brief Move constructor.
       */
      Base(Base&& other)
        : Copyable(std::move(other)),
          Identifiable(std::move(other))
      {}

      /**
       * @brief Destructor.
       */
      virtual ~Base() = default;

      /**
       * @brief Copy assignment is not allowed.
       */
      Base& operator=(const Base&) = delete;

      /**
       * @brief Move assignment is not allowed.
       */
      Base& operator=(Base&&) = delete;

      /**
       * @brief Gets the human-readable name of this object.
       * 
       * @return const char* Object name string, or nullptr if no name is set
       * 
       * Returns a string representation of the object type or expression.
       * Derived classes should override this method to provide meaningful
       * names for debugging and diagnostic purposes.
       */
      virtual Optional<StringView> getName() const
      {
        return {};
      }

      /**
       * @brief Materializes expression objects into concrete value objects.
       *
       * @tparam T Type of object to materialize
       * @param obj Object to materialize
       */
      template <class T>
      using DecayT = std::remove_cv_t<std::remove_reference_t<T>>;

      template <class T>
      static constexpr bool IsStackBackedObject =
        FormLanguage::RangeKindOf<DecayT<T>>::Value != FormLanguage::RangeKind::Unknown;

      /**
       * @brief Materializes stack-backed math objects.
       *
       * Supported categories are Boolean, Integer, Real, Complex,
       * SpatialVector and SpatialMatrix. Vector/matrix-like Eigen expressions
       * are explicitly materialized into spatial objects.
       */
      template <class T, std::enable_if_t<IsStackBackedObject<T>, int> = 0>
      constexpr
      auto object(T&& obj) const noexcept
      {
        using D = DecayT<T>;
        constexpr auto kind = FormLanguage::RangeKindOf<D>::Value;
        if constexpr (kind == FormLanguage::RangeKind::Boolean
                   || kind == FormLanguage::RangeKind::Integer
                   || kind == FormLanguage::RangeKind::Real
                   || kind == FormLanguage::RangeKind::Complex)
        {
          return static_cast<D>(std::forward<T>(obj));
        }
        else if constexpr (kind == FormLanguage::RangeKind::Vector)
        {
          using Scalar = typename D::Scalar;
          return Math::SpatialVector<Scalar>(std::forward<T>(obj));
        }
        else if constexpr (kind == FormLanguage::RangeKind::Matrix)
        {
          using Scalar = typename D::Scalar;
          return Math::SpatialMatrix<Scalar>(std::forward<T>(obj));
        }
        else
        {
          return static_cast<D>(std::forward<T>(obj));
        }
      }

      template <class T,
        std::enable_if_t<
          !IsStackBackedObject<T> && FormLanguage::IsPlainObject<std::remove_reference_t<T>>::Value,
        int> = 0>
      constexpr
      DecayT<T> object(T&& obj) const noexcept
      {
        return static_cast<DecayT<T>>(std::forward<T>(obj));
      }

      /**
       * @brief Forwards non-plain objects unchanged.
       * @tparam T Type of object (must not be a plain object type)
       * @param[in] obj Object to forward
       * @return Forwarded object
       *
       * This overload handles non-plain object types (such as scalars, references,
       * or expression templates) by forwarding them directly without storage.
       * It is selected via SFINAE when T is not a plain object.
       */
      template <class T,
        std::enable_if_t<
          !IsStackBackedObject<T> && !FormLanguage::IsPlainObject<std::remove_reference_t<T>>::Value,
        int> = 0>
      constexpr
      T object(T&& obj) const noexcept
      {
        return std::forward<T>(obj);
      }

      /**
       * @brief Clears internal temporary state.
       */
      void clear()
      {
        // No-op: object() now materializes values without heap storage.
      }

      /**
       * @brief Creates a polymorphic copy of this object.
       * @return Non-owning pointer to the copied object
       *
       * Pure virtual function that must be implemented by derived classes
       * to support polymorphic copying. The returned pointer is non-owning;
       * the caller is responsible for managing its lifetime.
       *
       * @note This is a CRTP function to be overridden in derived classes.
       */
      virtual Base* copy() const noexcept override = 0;

  };
}

#endif
