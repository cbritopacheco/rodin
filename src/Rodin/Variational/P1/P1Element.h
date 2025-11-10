/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_P1ELEMENT_H
#define RODIN_VARIATIONAL_P1_P1ELEMENT_H

/**
 * @file
 * @brief P1 (piecewise linear) finite element implementation.
 *
 * This file provides the P1Element class template for continuous piecewise
 * linear finite elements. P1 elements have:
 * - One DOF per vertex
 * - Linear basis functions: @f$ \phi_i(x_j) = \delta_{ij} @f$ (Lagrange property)
 * - Constant gradient per element: @f$ \nabla \phi_i|_K = \text{const} @f$
 *
 * P1 elements are the most common finite elements, providing first-order
 * convergence (@f$ O(h) @f$ for L² norm, @f$ O(h^2) @f$ for energy norm) and
 * forming the foundation for many FEM applications.
 */

#include <utility>

#include <boost/serialization/access.hpp>

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/FiniteElement.h"

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
   * @brief Continuous piecewise linear (degree 1) scalar Lagrange element.
   *
   * The P1Element provides a first-order finite element with:
   * - **DOF count**: One per vertex (@f$ n_v @f$ total)
   * - **Basis functions**: Linear functions satisfying @f$ \phi_i(x_j) = \delta_{ij} @f$
   * - **Gradient**: Constant on each element, @f$ \nabla \phi_i|_K = \text{const} @f$
   * - **Continuity**: Global C⁰ continuity across element interfaces
   *
   * P1 elements provide first-order convergence and are the standard choice
   * for elliptic PDEs like Poisson's equation @f$ -\Delta u = f @f$.
   *
   * @tparam Scalar Type of scalar range (e.g., Real, Complex)
   */
  template <class Scalar>
  class P1Element final : public FiniteElementBase<P1Element<Scalar>>
  {
    using G = Geometry::Polytope::Type;

    public:
      friend class boost::serialization::access;

      /// Parent class
      using Parent = FiniteElementBase<P1Element<Scalar>>;

      using ScalarType = Scalar;

      /// Type of range
      using RangeType = ScalarType;

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
          LinearForm(const LinearForm&) = delete;

          template <class T>
          ScalarType operator()(const T& v) const
          {
            const auto& vtx = Geometry::Polytope::Traits(m_g).getVertex(m_i);
            return v(vtx);
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
          using ReturnType = ScalarType;

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
              ReturnType operator()(const Math::SpatialPoint& r) const;

            private:
              const size_t m_i;
              const size_t m_local;
              const Geometry::Polytope::Type m_g;
          };

          class GradientFunction
          {
            public:
              using ReturnType = Math::SpatialVector<ScalarType>;

              constexpr
              GradientFunction(size_t local, Geometry::Polytope::Type g)
                : m_local(local), m_g(g)
              {}

              constexpr
              GradientFunction(const GradientFunction&) = default;

              const ReturnType& operator()(const Math::SpatialPoint& r) const
              {
                static thread_local ReturnType s_out;
                const size_t dim = Geometry::Polytope::Traits(m_g).getDimension();
                s_out.resize(dim);
                for (size_t i = 0; i < dim; ++i)
                  s_out(i) = DerivativeFunction<1>(i, m_local, m_g)(r);
                return s_out;
              }

            private:
              const size_t m_local;
              const Geometry::Polytope::Type m_g;
          };

          constexpr
          BasisFunction(size_t local, Geometry::Polytope::Type g)
            : m_local(local), m_g(g)
          {}

          constexpr
          BasisFunction(const BasisFunction&) = delete;

          constexpr
          ReturnType operator()(const Math::SpatialPoint& r) const;

          template <size_t Order>
          constexpr
          DerivativeFunction<Order> getDerivative(size_t i) const
          {
            return DerivativeFunction<Order>(i, m_local, m_g);
          }

          constexpr
          GradientFunction getGradient() const
          {
            return GradientFunction(m_local, m_g);
          }

        private:
          const size_t m_local;
          const Geometry::Polytope::Type m_g;
      };

      constexpr
      P1Element()
        : Parent(Geometry::Polytope::Type::Point)
      {}

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

      constexpr
      P1Element& operator=(const P1Element& other)
      {
        Parent::operator=(other);
        return *this;
      }

      constexpr
      P1Element& operator=(P1Element&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

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
      const Math::SpatialPoint& getNode(size_t i) const
      {
        return Geometry::Polytope::Traits(this->getGeometry()).getVertex(i);
      }

      const LinearForm& getLinearForm(size_t i) const
      {
        const Geometry::Polytope::Type g = this->getGeometry();
        switch (g)
        {
          case Geometry::Polytope::Type::Point:
          {
            static thread_local LinearForm s_lf(0, g);
            assert(i == 0);
            return s_lf;
          }
          case Geometry::Polytope::Type::Segment:
          {
            switch (i)
            {
              case 0:
              {
                static thread_local LinearForm s_lf(0, g);
                return s_lf;
              }
              case 1:
              {
                static thread_local LinearForm s_lf(1, g);
                return s_lf;
              }
              default:
              {
                assert(false);
                break;
              }
            }
            break;
          }
          case Geometry::Polytope::Type::Triangle:
          {
            switch (i)
            {
              case 0:
              {
                static thread_local LinearForm s_lf(0, g);
                return s_lf;
              }
              case 1:
              {
                static thread_local LinearForm s_lf(1, g);
                return s_lf;
              }
              case 2:
              {
                static thread_local LinearForm s_lf(2, g);
                return s_lf;
              }
              default:
              {
                assert(false);
                break;
              }
            }
            break;
          }
          case Geometry::Polytope::Type::Tetrahedron:
          {
            switch (i)
            {
              case 0:
              {
                static thread_local LinearForm s_lf(0, g);
                return s_lf;
              }
              case 1:
              {
                static thread_local LinearForm s_lf(1, g);
                return s_lf;
              }
              case 2:
              {
                static thread_local LinearForm s_lf(2, g);
                return s_lf;
              }
              case 3:
              {
                static thread_local LinearForm s_lf(3, g);
                return s_lf;
              }
              default:
              {
                assert(false);
                break;
              }
            }
            break;
          }
          case Geometry::Polytope::Type::Quadrilateral:
          {
            switch (i)
            {
              case 0:
              {
                static thread_local LinearForm s_lf(0, g);
                return s_lf;
              }
              case 1:
              {
                static thread_local LinearForm s_lf(1, g);
                return s_lf;
              }
              case 2:
              {
                static thread_local LinearForm s_lf(2, g);
                return s_lf;
              }
              case 3:
              {
                static thread_local LinearForm s_lf(3, g);
                return s_lf;
              }
              default:
              {
                assert(false);
                break;
              }
            }
            break;
          }
          case Geometry::Polytope::Type::Wedge:
          {
            switch (i)
            {
              case 0:
              {
                static thread_local LinearForm s_lf(0, g);
                return s_lf;
              }
              case 1:
              {
                static thread_local LinearForm s_lf(1, g);
                return s_lf;
              }
              case 2:
              {
                static thread_local LinearForm s_lf(2, g);
                return s_lf;
              }
              case 3:
              {
                static thread_local LinearForm s_lf(3, g);
                return s_lf;
              }
              case 4:
              {
                static thread_local LinearForm s_lf(4, g);
                return s_lf;
              }
              case 5:
              {
                static thread_local LinearForm s_lf(5, g);
                return s_lf;
              }
              default:
              {
                assert(false);
                break;
              }
            }
            break;
          }
        }
        static thread_local LinearForm s_null(0, g);
        assert(false);
        return s_null;
      }

      const BasisFunction& getBasis(size_t i) const
      {
        const Geometry::Polytope::Type g = this->getGeometry();
        switch (g)
        {
          case Geometry::Polytope::Type::Point:
          {
            static thread_local BasisFunction s_basis(0, g);
            assert(i == 0);
            return s_basis;
          }
          case Geometry::Polytope::Type::Segment:
          {
            switch (i)
            {
              case 0:
              {
                static thread_local BasisFunction s_basis(0, g);
                return s_basis;
              }
              case 1:
              {
                static thread_local BasisFunction s_basis(1, g);
                return s_basis;
              }
              default:
              {
                assert(false);
                break;
              }
            }
            break;
          }
          case Geometry::Polytope::Type::Triangle:
          {
            switch (i)
            {
              case 0:
              {
                static thread_local BasisFunction s_basis(0, g);
                return s_basis;
              }
              case 1:
              {
                static thread_local BasisFunction s_basis(1, g);
                return s_basis;
              }
              case 2:
              {
                static thread_local BasisFunction s_basis(2, g);
                return s_basis;
              }
              default:
              {
                assert(false);
                break;
              }
            }
            break;
          }
          case Geometry::Polytope::Type::Tetrahedron:
          {
            switch (i)
            {
              case 0:
              {
                static thread_local BasisFunction s_basis(0, g);
                return s_basis;
              }
              case 1:
              {
                static thread_local BasisFunction s_basis(1, g);
                return s_basis;
              }
              case 2:
              {
                static thread_local BasisFunction s_basis(2, g);
                return s_basis;
              }
              case 3:
              {
                static thread_local BasisFunction s_basis(3, g);
                return s_basis;
              }
              default:
              {
                assert(false);
                break;
              }
            }
            break;
          }
          case Geometry::Polytope::Type::Quadrilateral:
          {
            switch (i)
            {
              case 0:
              {
                static thread_local BasisFunction s_basis(0, g);
                return s_basis;
              }
              case 1:
              {
                static thread_local BasisFunction s_basis(1, g);
                return s_basis;
              }
              case 2:
              {
                static thread_local BasisFunction s_basis(2, g);
                return s_basis;
              }
              case 3:
              {
                static thread_local BasisFunction s_basis(3, g);
                return s_basis;
              }
              default:
              {
                assert(false);
                break;
              }
            }
            break;
          }
          case Geometry::Polytope::Type::Wedge:
          {
            switch (i)
            {
              case 0:
              {
                static thread_local BasisFunction s_basis(0, g);
                return s_basis;
              }
              case 1:
              {
                static thread_local BasisFunction s_basis(1, g);
                return s_basis;
              }
              case 2:
              {
                static thread_local BasisFunction s_basis(2, g);
                return s_basis;
              }
              case 3:
              {
                static thread_local BasisFunction s_basis(3, g);
                return s_basis;
              }
              case 4:
              {
                static thread_local BasisFunction s_basis(4, g);
                return s_basis;
              }
              case 5:
              {
                static thread_local BasisFunction s_basis(5, g);
                return s_basis;
              }
              default:
              {
                assert(false);
                break;
              }
            }
            break;
          }
        }
        static thread_local BasisFunction s_null(0, g);
        assert(false);
        return s_null;
      }

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

      template<class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        ar & boost::serialization::base_object<Parent>(*this);
      }
  };

  /**
   * @ingroup FiniteElements
   * @ingroup P1ElementSpecializations
   * @brief Continuous piecewise linear (degree 1) vector Lagrange element.
   *
   * Vector-valued P1 element with:
   * - **DOF count**: @f$ d \cdot n_v @f$ where @f$ d @f$ is vector dimension, @f$ n_v @f$ is vertex count
   * - **Basis functions**: @f$ \boldsymbol{\phi}_{i,j}(x) = \phi_i(x) \mathbf{e}_j @f$
   * - **Jacobian**: @f$ \mathbf{J}_{i,j} = \partial u_i/\partial x_j @f$ constant per element
   * - **Continuity**: C⁰ continuous vector field
   *
   * Used for elasticity, fluid mechanics, and vector-valued PDEs. Each component
   * uses P1 interpolation independently.
   *
   * @tparam Scalar Type of scalar components
   */
  template <class Scalar>
  class P1Element<Math::Vector<Scalar>> final
    : public FiniteElementBase<P1Element<Math::Vector<Scalar>>>
  {
    using G = Geometry::Polytope::Type;

    public:
      friend class boost::serialization::access;

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
          LinearForm(const LinearForm&) = delete;

          constexpr
          LinearForm(LinearForm&&) = default;

          template <class T>
          ScalarType operator()(const T& v) const
          {
            static thread_local RangeType s_out;
            const auto& vtx = P1Element<ScalarType>(m_g).getNode(m_local / m_vdim);
            s_out = v(vtx);
            return s_out.coeff(m_local % m_vdim);
          }

        private:
          const size_t m_vdim;
          const size_t m_local;
          const Geometry::Polytope::Type m_g;
      };

      class BasisFunction
      {
        public:
          using ReturnType = Math::Vector<ScalarType>;

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
              void operator()(Scalar& out, const Math::SpatialPoint& r) const
              {
                out = this->operator()(r);
              }

              constexpr
              Scalar operator()(const Math::SpatialPoint& rc) const
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

          class JacobianFunction
          {
            public:
              using ReturnType = Math::PointMatrix;

              constexpr
              JacobianFunction(size_t vdim, size_t local, Geometry::Polytope::Type g)
                : m_vdim(vdim), m_local(local), m_g(g)
              {}

              constexpr
              JacobianFunction(const JacobianFunction&) = default;

              constexpr
              JacobianFunction(JacobianFunction&&) = default;

              const ReturnType& operator()(const Math::SpatialPoint& r) const
              {
                static thread_local ReturnType s_out;
                const size_t dim = Geometry::Polytope::Traits(m_g).getDimension();
                s_out.resize(m_vdim, dim);
                for (size_t i = 0; i < m_vdim; ++i)
                {
                  for (size_t j = 0; j < dim; ++j)
                    s_out(i, j) = DerivativeFunction<1>(i, j, m_vdim, m_local, m_g)(r);
                }
                return s_out;
              }

            private:
              const size_t m_vdim;
              const size_t m_local;
              const Geometry::Polytope::Type m_g;
          };

          constexpr
          BasisFunction(size_t vdim, size_t local, Geometry::Polytope::Type g)
            : m_vdim(vdim), m_local(local), m_g(g)
          {}

          constexpr
          BasisFunction(const BasisFunction&) = delete;

          constexpr
          BasisFunction(BasisFunction&&) = default;

          const ReturnType& operator()(const Math::SpatialPoint& rc) const
          {
            static thread_local ReturnType s_out;
            s_out.resize(m_vdim);
            s_out.setZero();
            s_out.coeffRef(m_local % m_vdim) = P1Element<ScalarType>(m_g).getBasis(m_local / m_vdim)(rc);
            return s_out;
          }

          template <size_t Order>
          constexpr
          DerivativeFunction<Order> getDerivative(size_t i, size_t j) const
          {
            return DerivativeFunction<Order>(i, j, m_vdim, m_local, m_g);
          }

          constexpr
          JacobianFunction getJacobian() const
          {
            return JacobianFunction(m_vdim, m_local, m_g);
          }

        private:
          const size_t m_vdim;
          const size_t m_local;
          const Geometry::Polytope::Type m_g;
      };

      P1Element()
        : Parent(Geometry::Polytope::Type::Point), m_vdim(0)
      {}

      constexpr
      P1Element(size_t vdim, Geometry::Polytope::Type geometry)
        : Parent(geometry), m_vdim(vdim)
      {
        const size_t count = this->getCount();
        m_lfs.reserve(count);
        m_bs.reserve(count);
        for (size_t i = 0; i < count; ++i)
        {
          m_lfs.emplace_back(vdim, i, geometry);
          m_bs.emplace_back(vdim, i, geometry);
        }
      }

      constexpr
      P1Element(const P1Element& other)
        : Parent(other)
      {}

      constexpr
      P1Element(P1Element&& other)
        : Parent(std::move(other))
      {}

      constexpr
      P1Element& operator=(const P1Element& other)
      {
        Parent::operator=(other);
        m_vdim = other.m_vdim;
        m_lfs = other.m_lfs;
        m_bs = other.m_bs;
        return *this;
      }

      constexpr
      P1Element& operator=(P1Element&& other)
      {
        Parent::operator=(std::move(other));
        m_vdim = std::exchange(other.m_vdim, 0);
        m_lfs = std::move(other.m_lfs);
        m_bs = std::move(other.m_bs);
        return *this;
      }

      constexpr
      size_t getCount() const
      {
        return m_vdim * Geometry::Polytope::Traits(this->getGeometry()).getVertexCount();
      }

      constexpr
      const LinearForm& getLinearForm(size_t local) const
      {
        return m_lfs[local];
      }

      constexpr
      const BasisFunction& getBasis(size_t local) const
      {
        return m_bs[local];
      }

      constexpr
      const Math::SpatialPoint& getNode(size_t local) const
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

      template<class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        ar & boost::serialization::base_object<Parent>(*this);
        ar & m_vdim;
      }

    private:
      size_t m_vdim;

      std::vector<LinearForm> m_lfs;
      std::vector<BasisFunction> m_bs;
  };
}

#include "P1Element.hpp"

#endif
