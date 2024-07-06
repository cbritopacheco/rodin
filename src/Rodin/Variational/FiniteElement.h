/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENT_H
#define RODIN_VARIATIONAL_FINITEELEMENT_H

#include <unordered_map>

#include "Rodin/Math/Matrix.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/PolytopeTransformation.h"

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

      virtual ~FiniteElementBase() = default;

      constexpr
      FiniteElementBase& operator=(const FiniteElementBase& other)
      {
        m_g = other.m_g;
        return *this;
      }

      constexpr
      FiniteElementBase& operator=(FiniteElementBase&& other)
      {
        m_g = std::move(other.m_g);
        return *this;
      }

      inline
      constexpr
      Geometry::Polytope::Type getGeometry() const
      {
        return m_g;
      }

      /**
       * @brief Gets the number of degrees of freedom in the finite element.
       * @note CRTP method to be overriden in Derived class.
       */
      inline
      constexpr
      size_t getCount() const
      {
        return static_cast<const Derived&>(*this).getCount();
      }

      /**
       * @brief Gets the i-th degree of freedom on the finite element.
       */
      inline
      constexpr
      auto getNode(size_t i) const
      {
        return getNodes().col(i);
      }

      /**
       * @note CRTP method to be overriden in Derived class.
       */
      inline
      constexpr
      const Math::PointMatrix& getNodes() const
      {
        return static_cast<const Derived&>(*this).getNodes();
      }

      /**
       * @brief Gets the i-th basis function of the finite element.
       * @note CRTP method to be overriden in Derived class.
       */
      inline
      constexpr
      const auto& getBasis(size_t i) const
      {
        return static_cast<const Derived&>(*this).getBasis(i);
      }

      /**
       * @brief Gets the i-th linear function on the finite element.
       * @note CRTP method to be overriden in Derived class.
       */
      inline
      constexpr
      const auto& getLinearForm(size_t i) const
      {
        return static_cast<const Derived&>(*this).getLinearForm(i);
      }

      inline
      constexpr
      size_t getOrder() const
      {
        return static_cast<const Derived&>(*this).getOrder();
      }

    private:
      Geometry::Polytope::Type m_g;
  };
}

#endif

