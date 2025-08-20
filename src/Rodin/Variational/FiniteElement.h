/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENT_H
#define RODIN_VARIATIONAL_FINITEELEMENT_H

#include <boost/serialization/access.hpp>

#include "Rodin/Geometry/Polytope.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @defgroup FiniteElements Supported finite elements
   * @brief List of finite elements already implemented.
   */

  /**
   * @brief Base class for finite elements.
   */
  template <class Derived>
  class FiniteElementBase
  {
    public:
      using ScalarType = typename FormLanguage::Traits<Derived>::ScalarType;

      constexpr
      FiniteElementBase()
        : m_g(Geometry::Polytope::Type::Point)
      {}

      constexpr
      FiniteElementBase(Geometry::Polytope::Type g)
        : m_g(g)
      {}

      constexpr
      FiniteElementBase(const FiniteElementBase& other)
        : m_g(other.m_g)
      {}

      constexpr
      FiniteElementBase(FiniteElementBase&& other)
        : m_g(std::move(other.m_g))
      {}

      constexpr
      FiniteElementBase& operator=(const FiniteElementBase& other)
      {
        if (this != &other)
        {
          m_g = other.m_g;
        }
        return static_cast<Derived&>(*this);
      }

      constexpr
      FiniteElementBase& operator=(FiniteElementBase&& other)
      {
        if (this != &other)
        {
          m_g = std::move(other.m_g);
        }
        return static_cast<Derived&>(*this);
      }

      virtual ~FiniteElementBase() = default;

      constexpr
      Geometry::Polytope::Type getGeometry() const
      {
        return m_g;
      }

      /**
       * @brief Gets the number of degrees of freedom in the finite element.
       * @note CRTP method to be overriden in Derived class.
       */
      constexpr
      size_t getCount() const
      {
        return static_cast<const Derived&>(*this).getCount();
      }

      /**
       * @brief Gets the i-th degree of freedom on the finite element.
       */
      constexpr
      const Math::SpatialVector<Real>& getNode(size_t i) const
      {
        return static_cast<const Derived&>(*this).getNode(i);
      }

      /**
       * @brief Gets the i-th basis function of the finite element.
       * @note CRTP method to be overriden in Derived class.
       */
      constexpr
      decltype(auto) getBasis(size_t i) const
      {
        return static_cast<const Derived&>(*this).getBasis(i);
      }

      /**
       * @brief Gets the i-th linear function on the finite element.
       * @note CRTP method to be overriden in Derived class.
       */
      constexpr
      decltype(auto) getLinearForm(size_t i) const
      {
        return static_cast<const Derived&>(*this).getLinearForm(i);
      }

      constexpr
      size_t getOrder() const
      {
        return static_cast<const Derived&>(*this).getOrder();
      }

      template<class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        ar & m_g;
      }

    private:
      Geometry::Polytope::Type m_g;
  };
}

#endif

