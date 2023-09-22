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
 * @brief Header contatining definitions for the class P0Element.
 */

/**
 * @ingroup RodinDirectives
 * @brief Indicates the maximum vector dimension a P0Element
 */
#define RODIN_P0_MAX_VECTOR_DIMENSION 16

#include "Rodin/Types.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"
#include "Rodin/Geometry/GeometryIndexed.h"

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/FiniteElement.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"

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
   * @brief Degree 0 scalar Lagrange element
   */
  template <>
  class P0Element<Scalar> final : public FiniteElementBase<P0Element<Scalar>>
  {
    using G = Geometry::Polytope::Geometry;

    public:
      /// Type of range
      using RangeType = Scalar;

      /// Parent class
      using Parent = FiniteElementBase<P0Element<Scalar>>;

      /**
       * @brief Represents a linear form of a P0 scalar element.
       */
      class LinearForm
      {
        public:
          constexpr
          LinearForm(size_t i, Geometry::Polytope::Geometry g)
            : m_i(i), m_g(g)
          {
            assert(i < Geometry::Polytope::getVertexCount(g));
          }

          constexpr
          LinearForm(const LinearForm&) = default;

          constexpr
          LinearForm(LinearForm&&) = default;

          constexpr
          LinearForm& operator=(const LinearForm&) = default;

          constexpr
          LinearForm& operator=(LinearForm&&) = default;

          template <class T>
          inline
          constexpr
          auto operator()(const T& v) const
          {
            return v(s_nodes[m_g].col(m_i));
          }

        private:
          size_t m_i;
          Geometry::Polytope::Geometry m_g;
      };

      /**
       * @brief Represents a basis function of a P0 scalar element.
       */
      class BasisFunction
      {
        public:
          using ReturnType = Scalar;

          BasisFunction(size_t i, Geometry::Polytope::Geometry g)
            : m_i(i), m_g(g)
          {
            assert(i < Geometry::Polytope::getVertexCount(g));
          }

          BasisFunction(const BasisFunction&) = default;

          BasisFunction& operator=(const BasisFunction&) = default;

          BasisFunction& operator=(BasisFunction&&) = default;

          Scalar operator()(const Math::SpatialVector& r) const;

        private:
          size_t m_i;
          Geometry::Polytope::Geometry m_g;
      };

      /**
       * @brief Represents a gradient basis function of a P0 scalar element.
       */
      class GradientFunction
      {
        public:
          constexpr
          GradientFunction(size_t i, Geometry::Polytope::Geometry g)
            : m_i(i), m_g(g)
          {
            assert(i < Geometry::Polytope::getVertexCount(g));
          }

          constexpr
          GradientFunction(const GradientFunction&) = default;

          constexpr
          GradientFunction& operator=(const GradientFunction&) = default;

          constexpr
          GradientFunction& operator=(GradientFunction&&) = default;

          Math::SpatialVector operator()(const Math::SpatialVector& r) const
          {
            Math::SpatialVector out;
            operator()(out, r);
            return out;
          }

          void operator()(Math::SpatialVector& out, const Math::SpatialVector& r) const;

        private:
          size_t m_i;
          Geometry::Polytope::Geometry m_g;
      };

      constexpr
      P0Element() = default;

      constexpr
      P0Element(Geometry::Polytope::Geometry geometry)
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
      P0Element& operator=(const P0Element& other)
      {
        Parent::operator=(other);
        return *this;
      }

      constexpr
      P0Element& operator=(P0Element&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      /**
       * @brief Gets the number of degrees of freedom in the finite element.
       * @returns Number of degrees of freedom
       */
      inline
      size_t getCount() const
      {
        return Geometry::Polytope::getVertexCount(getGeometry());
      }

      inline
      constexpr
      const Math::Matrix& getNodes() const
      {
        return s_nodes[getGeometry()];
      }

      inline
      const auto& getLinearForm(size_t i) const
      {
        assert(i < getCount());
        return s_ls[getGeometry()][i];
      }

      inline
      const auto& getBasis(size_t i) const
      {
        assert(i < getCount());
        return s_basis[getGeometry()][i];
      }

      inline
      const auto& getGradient(size_t i) const
      {
        assert(i < getCount());
        return s_gradient[getGeometry()][i];
      }

    private:
      static const Geometry::GeometryIndexed<Math::Matrix> s_nodes;
      static const Geometry::GeometryIndexed<std::vector<LinearForm>> s_ls;
      static const Geometry::GeometryIndexed<std::vector<BasisFunction>> s_basis;
      static const Geometry::GeometryIndexed<std::vector<GradientFunction>> s_gradient;
  };
}

#endif

