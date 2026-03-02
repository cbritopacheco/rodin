/**
 * @file Base.h
 * @brief Base class for distance function computation models.
 *
 * This file provides the Base class template, which serves as a foundation
 * for various distance function computation methods by managing interior and
 * interface regions of the computational domain.
 */
#ifndef RODIN_MODELS_DITANCE_BASE_H
#define RODIN_MODELS_DITANCE_BASE_H

#include "Rodin/Geometry/Types.h"

namespace Rodin::Distance
{
  /**
   * @brief Base class for distance function computation models using CRTP.
   *
   * This class provides common functionality for distance function models,
   * specifically management of interior and interface regions. It uses the
   * Curiously Recurring Template Pattern (CRTP) to enable static polymorphism.
   *
   * @tparam Derived The derived class type (CRTP pattern)
   *
   * ## Usage
   * Derived classes should inherit from this base class to automatically
   * gain methods for setting and retrieving interior and interface regions:
   * ```cpp
   * class MyDistanceModel : public Base<MyDistanceModel>
   * {
   *   // Implementation
   * };
   * ```
   */
  template <class Derived>
  class Base
  {
    public:
      /**
       * @brief Sets the interior region attributes.
       *
       * @tparam A1 Type of first attribute
       * @tparam As Types of remaining attributes
       * @param[in] a1 First attribute to add to interior region
       * @param[in] as Additional attributes to add to interior region
       * @return Reference to the derived class for method chaining
       */
      template <class A1, class ... As>
      Derived& setInterior(A1&& a1, As&& ... as)
      {
        m_interior =
          FlatSet<Geometry::Attribute>{std::forward<A1>(a1), std::forward<As>(as)...};
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Sets the interior region attributes from a set.
       *
       * @param[in] interior Set of attributes defining the interior region
       * @return Reference to the derived class for method chaining
       */
      Derived& setInterior(const FlatSet<Geometry::Attribute>& interior)
      {
        m_interior = interior;
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Sets the interface region attributes.
       *
       * @tparam A1 Type of first attribute
       * @tparam As Types of remaining attributes
       * @param[in] a1 First attribute to add to interface region
       * @param[in] as Additional attributes to add to interface region
       * @return Reference to the derived class for method chaining
       */
      template <class A1, class ... As>
      Derived& setInterface(A1&& a1, As&& ... as)
      {
        m_interface =
          FlatSet<Geometry::Attribute>{std::forward<A1>(a1), std::forward<As>(as)...};
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Sets the interface region attributes from a set.
       *
       * @param[in] interface Set of attributes defining the interface region
       * @return Reference to the derived class for method chaining
       */
      Derived& setInterface(const FlatSet<Geometry::Attribute>& interface)
      {
        m_interface = interface;
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Gets the interior region attributes.
       *
       * @return Const reference to the set of interior region attributes
       */
      const auto& getInterior() const
      {
        return m_interior;
      }

      /**
       * @brief Gets the interface region attributes.
       *
       * @return Const reference to the set of interface region attributes
       */
      const auto& getInterface() const
      {
        return m_interface;
      }

    private:
      FlatSet<Geometry::Attribute> m_interior;   ///< Interior region attributes
      FlatSet<Geometry::Attribute> m_interface;  ///< Interface region attributes
  };
}

#endif

