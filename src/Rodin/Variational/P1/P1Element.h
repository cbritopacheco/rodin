/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_P1ELEMENT_H
#define RODIN_VARIATIONAL_P1_P1ELEMENT_H

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
  struct Traits<Variational::P1Element<Range>>
  {
    using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;
    using RangeType = Range;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup P1ElementSpecializations P1Element Template Specializations
   * @brief Template specializations of the P1Element class.
   * @see P1Element
   */

  /**
   * @ingroup FiniteElements
   * @ingroup P1ElementSpecializations
   * @brief Degree 1 scalar Lagrange element
   * @tparam Scalar Type of scalar range (e.g., Real, Complex)
   */
  template <class Scalar>
  class P1Element final : public FiniteElementBase<P1Element<Scalar>>
  {
    using G = Geometry::Polytope::Type;

    public:
      /// Parent class
      using Parent = FiniteElementBase<P1Element<Scalar>>;

      /// Type of range
      using RangeType = Scalar;

      /**
       * @brief Represents a linear form of a P1 scalar element.
       */
      class LinearForm
      {
        public:
          constexpr
          LinearForm(size_t i, Geometry::Polytope::Type g)
            : m_i(i), m_g(g)
          {}

          constexpr
          LinearForm(const LinearForm&) = default;

          template <class T>
          constexpr
          auto operator()(const T& v) const
          {
            return v(Geometry::Polytope::Traits(m_g).getVertex(m_i));
          }

        private:
          const size_t m_i;
          const Geometry::Polytope::Type m_g;
      };

      /**
       * @brief Represents a basis function of a P1 scalar element.
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
              DerivativeFunction(size_t i, size_t local, Geometry::Polytope::Type g)
                : m_i(i), m_local(local), m_g(g)
              {}

              constexpr
              DerivativeFunction(const DerivativeFunction&) = default;

              constexpr
              void operator()(ReturnType& out, const Math::SpatialVector<Real>& r) const
              {
                out = this->operator()(r);
              }

              constexpr
              ReturnType operator()(const Math::SpatialVector<Real>& r) const;

            private:
              const size_t m_i;
              const size_t m_local;
              const Geometry::Polytope::Type m_g;
          };

          constexpr
          BasisFunction(size_t local, Geometry::Polytope::Type g)
            : m_local(local), m_g(g)
          {}

          constexpr
          BasisFunction(const BasisFunction&) = default;

          constexpr
          void operator()(ReturnType& out, const Math::SpatialVector<Real>& r) const
          {
            out = this->operator()(r);
          }

          constexpr
          ReturnType operator()(const Math::SpatialVector<Real>& r) const;

          template <size_t Order>
          constexpr
          DerivativeFunction<Order> getDerivative(size_t i) const
          {
            return DerivativeFunction<Order>(i, m_local, m_g);
          }

        private:
          const size_t m_local;
          const Geometry::Polytope::Type m_g;
      };

      constexpr
      P1Element() = default;

      constexpr
      P1Element(Geometry::Polytope::Type geometry)
        : Parent(geometry)
      {}

      constexpr
      P1Element(const P1Element& other)
        : Parent(other)
      {}

      constexpr
      P1Element(P1Element&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Gets the number of degrees of freedom in the finite element.
       * @returns Number of degrees of freedom
       */
      constexpr
      size_t getCount() const
      {
        return Geometry::Polytope::Traits(this->getGeometry()).getVertexCount();
      }

      constexpr
      const Math::SpatialVector<Real>& getNode(size_t i) const
      {
        return Geometry::Polytope::Traits(this->getGeometry()).getVertex(i);
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
  };

  /**
   * @ingroup FiniteElements
   * @ingroup P1ElementSpecializations
   * @brief Degree 1 vector Lagrange element
   */
  template <class Scalar>
  class P1Element<Math::Vector<Scalar>> final
    : public FiniteElementBase<P1Element<Math::Vector<Scalar>>>
  {
    using G = Geometry::Polytope::Type;

    public:
      /// Parent class
      using Parent = FiniteElementBase<P1Element<Math::Vector<Scalar>>>;

      using ScalarType = Scalar;

      /// Type of range
      using RangeType = Math::Vector<Scalar>;

      class LinearForm
      {
        public:
          constexpr
          LinearForm(size_t vdim, size_t local, Geometry::Polytope::Type g)
            : m_vdim(vdim), m_local(local), m_g(g)
          {}

          constexpr
          LinearForm(const LinearForm&) = default;

          constexpr
          LinearForm(LinearForm&&) = default;

          template <class T>
          constexpr
          auto operator()(const T& v) const
          {
            return v(P1Element<ScalarType>(m_g).getNode(m_local / m_vdim)).coeff(m_local % m_vdim);
          }

        private:
          const size_t m_vdim;
          const size_t m_local;
          const Geometry::Polytope::Type m_g;
      };

      class BasisFunction
      {
        public:
          /**
           * @brief Represents a derivative function of a P1 vector element.
           * @tparam Order Order of the derivative (0 for function, 1 for first
           * derivative, etc.)
           */
          template <size_t Order>
          class DerivativeFunction
          {
            public:
              /**
               * @brief Constructs a derivative function for the P1 element.
               * @param i Index of the vector component
               * @param j Index of the space variable
               * @param local Local index of the basis function
               * @param g Geometry type of the polytope
               */
              constexpr
              DerivativeFunction(size_t i, size_t j, size_t vdim, size_t local, Geometry::Polytope::Type g)
                : m_i(i), m_j(j), m_vdim(vdim), m_local(local), m_g(g)
              {}

              constexpr
              DerivativeFunction(const DerivativeFunction&) = default;

              constexpr
              void operator()(Scalar& out, const Math::SpatialVector<Real>& r) const
              {
                out = this->operator()(r);
              }

              constexpr
              Scalar operator()(const Math::SpatialVector<Real>& rc) const
              {
                if constexpr (Order == 0)
                {
                  if (m_i == m_local % m_vdim)
                  {
                    return P1Element<ScalarType>(m_g).getBasis(m_local / m_vdim)(rc);
                  }
                  else
                  {
                    return 0;
                  }
                }
                else if constexpr (Order == 1)
                {
                  if (m_i == m_local % m_vdim)
                  {
                    return P1Element<ScalarType>(m_g).getBasis(m_local / m_vdim).template getDerivative<1>(m_j)(rc);
                  }
                  else
                  {
                    return 0;
                  }
                }
                else
                {
                  return 0;
                }
              }

            private:
              const size_t m_i, m_j;
              const size_t m_vdim;
              const size_t m_local;
              const Geometry::Polytope::Type m_g;
          };

          constexpr
          BasisFunction(size_t vdim, size_t local, Geometry::Polytope::Type g)
            : m_vdim(vdim), m_local(local), m_g(g)
          {}

          constexpr
          BasisFunction(const BasisFunction&) = default;

          constexpr
          BasisFunction(BasisFunction&&) = default;

          Math::Vector<ScalarType> operator()(const Math::SpatialVector<Real>& r) const
          {
            Math::Vector<ScalarType> res;
            operator()(res, r);
            return res;
          }

          constexpr
          void operator()(Math::Vector<ScalarType>& out, const Math::SpatialVector<Real>& rc) const
          {
            out.resize(m_vdim);
            out.setZero();
            out.coeffRef(m_local % m_vdim) = P1Element<ScalarType>(m_g).getBasis(m_local / m_vdim)(rc);
          }

          template <size_t Order>
          constexpr
          DerivativeFunction<Order> getDerivative(size_t i, size_t j) const
          {
            return DerivativeFunction<Order>(i, j, m_vdim, m_local, m_g);
          }

        private:
          const size_t m_vdim;
          const size_t m_local;
          const Geometry::Polytope::Type m_g;
      };

      P1Element() = default;

      constexpr
      P1Element(size_t vdim, Geometry::Polytope::Type geometry)
        : Parent(geometry), m_vdim(vdim)
      {}

      constexpr
      P1Element(const P1Element& other)
        : Parent(other)
      {}

      constexpr
      P1Element(P1Element&& other)
        : Parent(std::move(other))
      {}

      constexpr
      size_t getCount() const
      {
        return m_vdim * Geometry::Polytope::Traits(this->getGeometry()).getVertexCount();
      }

      constexpr
      LinearForm getLinearForm(size_t local) const
      {
        return LinearForm(m_vdim, local, this->getGeometry());
      }

      constexpr
      BasisFunction getBasis(size_t local) const
      {
        return BasisFunction(m_vdim, local, this->getGeometry());
      }

      constexpr
      const Math::SpatialVector<Real>& getNode(size_t local) const
      {
        return P1Element<ScalarType>(this->getGeometry()).getNode(local / m_vdim);
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
      const size_t m_vdim;
  };
}

#include "P1Element.hpp"

#endif
