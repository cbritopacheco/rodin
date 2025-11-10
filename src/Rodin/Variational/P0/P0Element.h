/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0_P0ELEMENT_H
#define RODIN_VARIATIONAL_P0_P0ELEMENT_H

/**
 * @file
 * @brief P0 (piecewise constant) finite element implementation.
 *
 * This file provides the P0Element class template for piecewise constant
 * finite elements. P0 elements have:
 * - One degree of freedom per element (at the barycenter)
 * - Constant basis function: @f$ \phi(x) = 1 @f$
 * - Zero gradient: @f$ \nabla \phi = 0 @f$
 *
 * P0 elements are commonly used in discontinuous Galerkin (DG) methods,
 * mixed finite element formulations, and as test spaces for finite volume
 * schemes.
 */

#include "Rodin/Types.h"

#include "Rodin/Math/Traits.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"
#include "Rodin/Geometry/GeometryIndexed.h"

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/FiniteElement.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"

namespace Rodin::FormLanguage
{
  /**
   * @ingroup TraitsSpecializations
   */
  template <class Range>
  struct Traits<Variational::P0Element<Range>>
  {
    using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;
    using RangeType = Range;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup P0ElementSpecializations P0Element Template Specializations
   * @brief Template specializations of the P0Element class.
   * @see P0Element
   */

  /**
   * @ingroup FiniteElements
   * @ingroup P0ElementSpecializations
   * @brief Piecewise constant (degree 0) scalar Lagrange element.
   *
   * The P0Element represents a piecewise constant finite element with:
   * - **DOF count**: 1 per element (located at barycenter)
   * - **Basis function**: @f$ \phi(x) = 1 @f$ for all @f$ x @f$ in element
   * - **Derivatives**: @f$ \nabla \phi = 0 @f$ (constant function has zero gradient)
   *
   * P0 elements are discontinuous across element interfaces, making them
   * suitable for DG methods, flux computations, and element-wise constant
   * approximations.
   *
   * @tparam Scalar Type of scalar range (e.g., Real, Complex)
   */
  template <class Scalar>
  class P0Element final : public FiniteElementBase<P0Element<Scalar>>
  {
    using G = Geometry::Polytope::Type;

    public:
      /// Parent class
      using Parent = FiniteElementBase<P0Element<Scalar>>;

      using ScalarType = Scalar;

      /// Type of range
      using RangeType = Scalar;

      /**
       * @brief Represents a linear form of a P0 scalar element.
       */
      class LinearForm
      {
        public:
          constexpr
          LinearForm(Geometry::Polytope::Type g)
            : m_g(g)
          {}

          constexpr
          LinearForm(const LinearForm&) = default;

          template <class T>
          constexpr
          ScalarType operator()(const T& v) const
          {
            return v(s_nodes[m_g]);
          }

        private:
          const Geometry::Polytope::Type m_g;
      };

      /**
       * @brief Represents a basis function of a P0 scalar element.
       */
      class BasisFunction
      {
        public:
          using ReturnType = Scalar;

          template <size_t Order>
          class DerivativeFunction
          {
            public:
              constexpr
              DerivativeFunction() = default;

              constexpr
              DerivativeFunction(const DerivativeFunction&) = default;

              constexpr
              ReturnType operator()(const Math::SpatialVector<Real>& r) const
              {
                return 0;
              }
          };

          constexpr
          BasisFunction() = default;

          constexpr
          BasisFunction(const BasisFunction&) = default;

          constexpr
          ReturnType operator()(const Math::SpatialVector<Real>& r) const
          {
            return 1;
          }

          template <size_t Order>
          constexpr
          DerivativeFunction<Order> getDerivative(size_t) const
          {
            return DerivativeFunction<Order>();
          }
      };

      constexpr
      P0Element() = default;

      constexpr
      P0Element(Geometry::Polytope::Type geometry)
        : Parent(geometry)
      {}

      constexpr
      P0Element(const P0Element& other)
        : Parent(other)
      {}

      constexpr
      P0Element(P0Element&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Gets the number of degrees of freedom in the finite element.
       * @returns Number of degrees of freedom
       */
      constexpr
      size_t getCount() const
      {
        return 1;
      }

      constexpr
      const Math::SpatialVector<Real>& getNode(size_t i) const
      {
        return s_nodes[this->getGeometry()][i];
      }

      constexpr
      LinearForm getLinearForm(size_t) const
      {
        return LinearForm(this->getGeometry());
      }

      constexpr
      BasisFunction getBasis(size_t) const
      {
        return BasisFunction();
      }

      constexpr
      size_t getOrder() const
      {
        return 0;
      }

    private:
      static const Geometry::GeometryIndexed<Math::SpatialVector<Real>> s_nodes;
  };

  /**
   * @ingroup FiniteElements
   * @ingroup P0ElementSpecializations
   * @brief Piecewise constant (degree 0) vector Lagrange element.
   *
   * Vector-valued P0 element with:
   * - **DOF count**: @f$ d @f$ per element (one for each vector component)
   * - **Basis functions**: @f$ \boldsymbol{\phi}_i(x) = \mathbf{e}_i @f$ (unit vectors)
   * - **Derivatives**: @f$ \nabla \boldsymbol{\phi}_i = \mathbf{0} @f$ (zero gradient)
   *
   * Used for vector-valued DG approximations (e.g., velocity in compressible flow).
   *
   * @tparam Scalar Type of scalar components
   */
  template <class Scalar>
  class P0Element<Math::Vector<Scalar>> final
    : public FiniteElementBase<P0Element<Math::Vector<Scalar>>>
  {
    using G = Geometry::Polytope::Type;

    public:
      /// Parent class
      using Parent = FiniteElementBase<P0Element>;

      using ScalarType = Scalar;

      /// Type of range
      using RangeType = Math::Vector<Scalar>;

      class LinearForm
      {
        public:
          constexpr
          LinearForm()
            : m_local(0), m_g(Geometry::Polytope::Type::Point)
          {}

          constexpr
          LinearForm(size_t local, Geometry::Polytope::Type g)
            : m_local(local), m_g(g)
          {}

          constexpr
          LinearForm(const LinearForm&) = default;

          constexpr
          LinearForm(LinearForm&&) = default;

          template <class T>
          constexpr
          decltype(auto) operator()(const T& v) const
          {
            const size_t vdim = Geometry::Polytope::Traits(m_g).getDimension();
            return v(P0Element<ScalarType>(m_g).getNode(m_local / vdim)).coeff(m_local % vdim);
          }

        private:
          const size_t m_local;
          const Geometry::Polytope::Type m_g;
      };

      class BasisFunction
      {
        public:
          template <size_t Order>
          class DerivativeFunction
          {
            public:
              constexpr
              DerivativeFunction() = default;

              constexpr
              DerivativeFunction(const DerivativeFunction&) = default;

              constexpr
              ScalarType operator()(const Math::SpatialVector<Real>& rc) const
              {
                return ScalarType(0);
              }
          };

          constexpr
          BasisFunction()
            : m_local(0), m_g(Geometry::Polytope::Type::Point)
          {}

          constexpr
          BasisFunction(size_t local, Geometry::Polytope::Type g)
            : m_local(local), m_g(g)
          {}

          constexpr
          BasisFunction(const BasisFunction&) = default;

          constexpr
          BasisFunction(BasisFunction&&) = default;

          auto operator()(const Math::SpatialVector<ScalarType>& r) const
          {
            const size_t vdim = Geometry::Polytope::Traits(m_g).getDimension();
            return Math::Vector<ScalarType>::Zero(vdim);
          }

          template <size_t Order>
          constexpr
          DerivativeFunction<Order> getDerivative(size_t i, size_t j) const
          {
            return DerivativeFunction<Order>(i, m_local, m_g);
          }

        private:
          const size_t m_local;
          const Geometry::Polytope::Type m_g;
      };

      P0Element() = default;

      constexpr
      P0Element(Geometry::Polytope::Type geometry)
        : Parent(geometry)
      {}

      constexpr
      P0Element(const P0Element& other)
        : Parent(other)
      {}

      constexpr
      P0Element(P0Element&& other)
        : Parent(std::move(other))
      {}

      constexpr
      size_t getCount() const
      {
        return Geometry::Polytope::Traits(this->getGeometry()).getDimension();
      }

      constexpr
      auto getLinearForm(size_t local) const
      {
        return LinearForm(local, this->getGeometry());
      }

      constexpr
      BasisFunction getBasis(size_t local) const
      {
        return BasisFunction(local, this->getGeometry());
      }

      constexpr
      const Math::SpatialVector<Real>& getNode(size_t local) const
      {
        const size_t vdim = Geometry::Polytope::Traits(this->getGeometry()).getDimension();
        return P0Element<ScalarType>(this->getGeometry()).getNode(local / vdim);
      }

      constexpr
      size_t getOrder() const
      {
        return 0;
      }
  };
}

#include "P0Element.hpp"

#endif
