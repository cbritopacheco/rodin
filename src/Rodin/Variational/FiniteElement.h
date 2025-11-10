/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file FiniteElement.h
 * @brief Base classes for finite element definitions.
 *
 * This file defines the FiniteElementBase class, which provides the foundation
 * for all finite element types. A finite element is characterized by the triple
 * @f$ (K, P, \Sigma) @f$ defining its geometry, polynomial space, and DOFs.
 *
 * ## Mathematical Foundation
 * A finite element is defined by:
 * - @f$ K @f$: Reference element geometry (interval, triangle, tetrahedron, etc.)
 * - @f$ P @f$: Polynomial space (e.g., @f$ \mathbb{P}_k @f$, @f$ \mathbb{Q}_k @f$)
 * - @f$ \Sigma @f$: Degrees of freedom (nodal values, moments, etc.)
 *
 * ## Element Types
 * Common finite element families:
 * - **Lagrange**: @f$ C^0 @f$ continuous, nodal DOFs at vertices/edges/faces
 * - **Hermite**: @f$ C^1 @f$ continuous, includes derivative DOFs
 * - **Raviart-Thomas**: @f$ H(\text{div}) @f$ conforming
 * - **Nédélec**: @f$ H(\text{curl}) @f$ conforming
 * - **Discontinuous**: @f$ L^2 @f$ elements for DG methods
 *
 * ## Polynomial Spaces
 * - @f$ \mathbb{P}_k @f$: Complete polynomials of degree @f$ \leq k @f$
 * - @f$ \mathbb{Q}_k @f$: Tensor product polynomials of degree @f$ \leq k @f$ per direction
 * - Serendipity: Subset of @f$ \mathbb{Q}_k @f$ with reduced DOFs
 *
 * ## Reference Elements
 * Standard geometries:
 * - 1D: Interval [0,1]
 * - 2D: Triangle (vertices at (0,0), (1,0), (0,1)) or quadrilateral [0,1]²
 * - 3D: Tetrahedron, hexahedron, prism, pyramid
 *
 * @see P0, P1, FiniteElementSpace
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

