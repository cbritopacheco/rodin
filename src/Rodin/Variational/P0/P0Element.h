/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0_P0ELEMENT_H
#define RODIN_VARIATIONAL_P0_P0ELEMENT_H

#include "Rodin/Types.h"

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
   * @brief Degree 1 scalar Lagrange element
   * @see @m_defelement{Lagrange,https://defelement.com/elements/lagrange.html}
   */
  template <class Scalar>
  class P0Element final : public FiniteElementBase<P0Element<Scalar>>
  {
    using G = Geometry::Polytope::Type;

    public:
      /// Parent class
      using Parent = FiniteElementBase<P0Element<Scalar>>;

      /// Type of range
      using RangeType = Scalar;

      /**
       * @brief Represents a linear form of a P0 scalar element.
       */
      class LinearForm
      {
        public:
          constexpr
          LinearForm(size_t i, Geometry::Polytope::Type g)
            : m_i(i), m_g(g)
          {
            assert(i < Geometry::Polytope::getVertexCount(g));
          }

          constexpr
          LinearForm(const LinearForm&) = default;

          template <class T>
          constexpr
          auto operator()(const T& v) const
          {
            return v(s_nodes[m_g][m_i]);
          }

        private:
          const size_t m_i;
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
              void operator()(ReturnType& out, const Math::SpatialVector<Real>& r) const
              {
                out = 0;
              }

              constexpr
              ReturnType operator()(const Math::SpatialVector<Real>& r) const
              {
                return 0;
              }
          };

          constexpr
          BasisFunction(size_t local, Geometry::Polytope::Type g)
            : m_local(local), m_g(g)
          {
            assert(local < Geometry::Polytope::getVertexCount(g));
          }

          constexpr
          BasisFunction(const BasisFunction&) = default;

          constexpr
          void operator()(ReturnType& out, const Math::SpatialVector<Real>& r) const
          {
            out = 1;
          }

          constexpr
          ReturnType operator()(const Math::SpatialVector<Real>& r) const
          {
            return 1;
          }

          template <size_t Order>
          constexpr
          DerivativeFunction<Order> getDerivative(size_t i) const
          {
            return DerivativeFunction<Order>();
          }

        private:
          const size_t m_local;
          const Geometry::Polytope::Type m_g;
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
      LinearForm getLinearForm(size_t i) const
      {
        return LinearForm(i, this->getGeometry());
      }

      constexpr
      BasisFunction getBasis(size_t i) const
      {
        return BasisFunction(i, this->getGeometry());
      }

      constexpr
      size_t getOrder() const
      {
        switch (this->getGeometry())
        {
          case Geometry::Polytope::Type::Point:
            return 0;
          case Geometry::Polytope::Type::Segment:
          case Geometry::Polytope::Type::Triangle:
          case Geometry::Polytope::Type::Tetrahedron:
            return 1;
          case Geometry::Polytope::Type::Quadrilateral:
          case Geometry::Polytope::Type::Wedge:
            return 2;
        }
        assert(false);
        return 0;
      }

    private:
      static const Geometry::GeometryIndexed<std::vector<Math::SpatialVector<Real>>> s_nodes;
  };

  /**
   * @ingroup FiniteElements
   * @ingroup P0ElementSpecializations
   * @brief Degree 1 vector Lagrange element
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
          auto operator()(const T& v) const
          {
            const size_t vdim = Geometry::Polytope::getGeometryDimension(m_g);
            const P0Element<ScalarType> P0e(m_g);
            return v(P0e.getNode(m_local / vdim)).coeff(m_local % vdim);
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
              void operator()(Scalar& out, const Math::SpatialVector<Real>& r) const
              {
                out = 0;
              }

              constexpr
              Scalar operator()(const Math::SpatialVector<Real>& rc) const
              {
                return 0;
              }

            private:
              const size_t m_i, m_j;
              const size_t m_local;
              const Geometry::Polytope::Type m_g;
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

          Math::Vector<ScalarType> operator()(const Math::SpatialVector<ScalarType>& r) const
          {
            Math::Vector<ScalarType> res;
            operator()(res, r);
            return res;
          }

          constexpr
          void operator()(Math::Vector<ScalarType>& out, const Math::SpatialVector<Real>& rc) const
          {
            const size_t vdim = Geometry::Polytope::getGeometryDimension(m_g);
            out.resize(vdim);
            out.setZero();
            out.coeffRef(m_local % vdim) = P0Element<ScalarType>(m_g).getBasis(m_local / vdim)(rc);
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
        return Geometry::Polytope::getVertexCount(this->getGeometry()
            ) * Geometry::Polytope::getGeometryDimension(this->getGeometry());
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
        const size_t vdim = Geometry::Polytope::getGeometryDimension(this->getGeometry());
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
