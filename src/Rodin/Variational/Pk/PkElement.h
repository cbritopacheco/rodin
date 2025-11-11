/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PK_PKELEMENT_H
#define RODIN_VARIATIONAL_PK_PKELEMENT_H

/**
 * @file
 * @brief Pk (piecewise polynomial degree k) finite element implementation.
 *
 * This file provides the PkElement class template for continuous piecewise
 * polynomial finite elements of arbitrary degree k. Pk elements have:
 * - Lagrange basis functions of degree k
 * - DOFs at Lagrange nodes (vertices, edge points, face points, volume points)
 * - Polynomial gradient of degree k-1
 *
 * Pk elements provide k-th order convergence and generalize P0 (k=0) and
 * P1 (k=1) elements.
 */

#include <cstddef>
#include <array>
#include <vector>
#include <utility>

#include <boost/serialization/access.hpp>

#include "Rodin/Types.h"
#include "Rodin/Math/Traits.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"
#include "Rodin/Geometry/Polytope.h"

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/FiniteElement.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"

/**
 * @ingroup RodinDirectives
 * @brief Indicates the maximum vector dimension a PkElement
 */
#define RODIN_PK_MAX_VECTOR_DIMENSION 16

namespace Rodin::FormLanguage
{
  /**
   * @ingroup TraitsSpecializations
   */
  template <size_t K, class Range>
  struct Traits<Variational::PkElement<K, Range>>
  {
    using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;
    using RangeType = Range;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup PkElementSpecializations PkElement Template Specializations
   * @brief Template specializations of the PkElement class.
   * @see PkElement
   */

  /**
   * @ingroup FiniteElements
   * @ingroup PkElementSpecializations
   * @brief Continuous piecewise polynomial (degree k) scalar Lagrange element.
   *
   * The PkElement provides a k-th order finite element with:
   * - **DOF count**: Depends on geometry and degree k
   * - **Basis functions**: Lagrange polynomials of degree k satisfying @f$ \phi_i(x_j) = \delta_{ij} @f$
   * - **Gradient**: Polynomial of degree k-1
   * - **Continuity**: Global C⁰ continuity across element interfaces
   *
   * Pk elements provide k-th order convergence for smooth solutions.
   *
   * @tparam K Polynomial degree
   * @tparam Scalar Type of scalar range (e.g., Real, Complex)
   */
  template <size_t K, class Scalar>
  class PkElement final : public FiniteElementBase<PkElement<K, Scalar>>
  {
    using G = Geometry::Polytope::Type;

    public:
      friend class boost::serialization::access;

      /// Parent class
      using Parent = FiniteElementBase<PkElement<K, Scalar>>;

      using ScalarType = Scalar;

      /// Type of range
      using RangeType = ScalarType;

      /**
       * @brief Represents a linear form of a Pk scalar element.
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
          ScalarType operator()(const T& v) const
          {
            const auto& node = PkElement<K, Scalar>(m_g).getNode(m_i);
            return v(node);
          }

        private:
          const size_t m_i;
          const Geometry::Polytope::Type m_g;
      };

      /**
       * @brief Represents a basis function of a Pk scalar element.
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
          BasisFunction(const BasisFunction&) = default;

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

      PkElement()
        : Parent(Geometry::Polytope::Type::Point)
      {}

      PkElement(Geometry::Polytope::Type geometry)
        : Parent(geometry)
      {
        buildNodes();
      }

      constexpr
      PkElement(const PkElement& other)
        : Parent(other), m_nodes(other.m_nodes)
      {}

      constexpr
      PkElement(PkElement&& other)
        : Parent(std::move(other)), m_nodes(std::move(other.m_nodes))
      {}

      constexpr
      PkElement& operator=(const PkElement& other)
      {
        Parent::operator=(other);
        m_nodes = other.m_nodes;
        return *this;
      }

      constexpr
      PkElement& operator=(PkElement&& other)
      {
        Parent::operator=(std::move(other));
        m_nodes = std::move(other.m_nodes);
        return *this;
      }

      /**
       * @brief Gets the number of degrees of freedom in the finite element.
       * @returns Number of degrees of freedom
       */
      constexpr
      size_t getCount() const;

      constexpr
      const Math::SpatialPoint& getNode(size_t i) const
      {
        return m_nodes[i];
      }

      const LinearForm& getLinearForm(size_t i) const;

      const BasisFunction& getBasis(size_t i) const;

      constexpr
      size_t getOrder() const
      {
        return K;
      }

      template<class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        ar & boost::serialization::base_object<Parent>(*this);
      }

    private:
      void buildNodes();

      std::vector<Math::SpatialPoint> m_nodes;
  };

  /**
   * @ingroup FiniteElements
   * @ingroup PkElementSpecializations
   * @brief Continuous piecewise polynomial (degree k) vector Lagrange element.
   *
   * Vector-valued Pk element with:
   * - **DOF count**: @f$ d \times @f$ (number of scalar DOFs) where @f$ d @f$ is vector dimension
   * - **Basis functions**: @f$ \boldsymbol{\phi}_{i,j}(x) = \phi_i(x) \mathbf{e}_j @f$
   * - **Jacobian**: @f$ \mathbf{J}_{i,j} = \partial u_i/\partial x_j @f$
   * - **Continuity**: C⁰ continuous vector field
   *
   * Used for elasticity, fluid mechanics, and vector-valued PDEs. Each component
   * uses Pk interpolation independently.
   *
   * @tparam K Polynomial degree
   * @tparam Scalar Type of scalar components
   */
  template <size_t K, class Scalar>
  class PkElement<K, Math::Vector<Scalar>> final
    : public FiniteElementBase<PkElement<K, Math::Vector<Scalar>>>
  {
    using G = Geometry::Polytope::Type;

    public:
      friend class boost::serialization::access;

      /// Parent class
      using Parent = FiniteElementBase<PkElement<K, Math::Vector<Scalar>>>;

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
          ScalarType operator()(const T& v) const
          {
            static thread_local RangeType s_out;
            const auto& node = PkElement<K, ScalarType>(m_g).getNode(m_local / m_vdim);
            s_out = v(node);
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
           * @brief Represents a derivative function of a Pk vector element.
           * @tparam Order Order of the derivative (0 for function, 1 for first
           * derivative, etc.)
           */
          template <size_t Order>
          class DerivativeFunction
          {
            public:
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
                    return PkElement<K, ScalarType>(m_g).getBasis(m_local / m_vdim)(rc);
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
                    return PkElement<K, ScalarType>(m_g).getBasis(m_local / m_vdim).template getDerivative<1>(m_j)(rc);
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
          BasisFunction(const BasisFunction&) = default;

          constexpr
          BasisFunction(BasisFunction&&) = default;

          const ReturnType& operator()(const Math::SpatialPoint& rc) const
          {
            static thread_local ReturnType s_out;
            s_out.resize(m_vdim);
            s_out.setZero();
            s_out.coeffRef(m_local % m_vdim) = PkElement<K, ScalarType>(m_g).getBasis(m_local / m_vdim)(rc);
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

      PkElement()
        : Parent(Geometry::Polytope::Type::Point), m_vdim(0)
      {}

      constexpr
      PkElement(size_t vdim, Geometry::Polytope::Type geometry)
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
      PkElement(const PkElement& other)
        : Parent(other), m_vdim(other.m_vdim), m_lfs(other.m_lfs), m_bs(other.m_bs)
      {}

      constexpr
      PkElement(PkElement&& other)
        : Parent(std::move(other)), m_vdim(std::move(other.m_vdim)),
          m_lfs(std::move(other.m_lfs)), m_bs(std::move(other.m_bs))
      {}

      constexpr
      PkElement& operator=(const PkElement& other)
      {
        Parent::operator=(other);
        m_vdim = other.m_vdim;
        m_lfs = other.m_lfs;
        m_bs = other.m_bs;
        return *this;
      }

      constexpr
      PkElement& operator=(PkElement&& other)
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
        return m_vdim * PkElement<K, ScalarType>(this->getGeometry()).getCount();
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
        return PkElement<K, ScalarType>(this->getGeometry()).getNode(local / m_vdim);
      }

      constexpr
      size_t getOrder() const
      {
        return K;
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

#include "PkElement.hpp"

#endif

