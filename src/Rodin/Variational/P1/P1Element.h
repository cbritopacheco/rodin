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
 * @brief Header contatining definitions for the class P1Element.
 */

/**
 * @ingroup RodinDirectives
 * @brief Indicates the maximum vector dimension a P1Element
 */
#define RODIN_P1_MAX_VECTOR_DIMENSION 16

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
  struct Traits<Variational::P1Element<Range>>
  {
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
   * @see @m_defelement{Lagrange,https://defelement.com/elements/lagrange.html}
   */
  template <>
  class P1Element<Scalar> final : public FiniteElementBase<P1Element<Scalar>>
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
          Geometry::Polytope::Type m_g;
      };

      /**
       * @brief Represents a basis function of a P1 scalar element.
       */
      class BasisFunction
      {
        public:
          using ReturnType = Scalar;

          BasisFunction(size_t i, Geometry::Polytope::Type g)
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
          Geometry::Polytope::Type m_g;
      };

      /**
       * @brief Represents a gradient basis function of a P1 scalar element.
       */
      class GradientFunction
      {
        public:
          constexpr
          GradientFunction(size_t i, Geometry::Polytope::Type g)
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
          Geometry::Polytope::Type m_g;
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

      inline
      size_t getOrder() const
      {
        switch (getGeometry())
        {
          case Geometry::Polytope::Type::Point:
            return 0;
          case Geometry::Polytope::Type::Segment:
          case Geometry::Polytope::Type::Triangle:
          case Geometry::Polytope::Type::Tetrahedron:
            return 1;
          case Geometry::Polytope::Type::Quadrilateral:
            return 2;
        }
        assert(false);
        return 0;
      }

    private:
      static const Geometry::GeometryIndexed<Math::Matrix> s_nodes;
      static const Geometry::GeometryIndexed<std::vector<LinearForm>> s_ls;
      static const Geometry::GeometryIndexed<std::vector<BasisFunction>> s_basis;
      static const Geometry::GeometryIndexed<std::vector<GradientFunction>> s_gradient;
  };

  /**
   * @ingroup FiniteElements
   * @ingroup P1ElementSpecializations
   * @brief Degree 1 vector Lagrange element
   *
   * @see @m_defelement{Vector Lagrange,https://defelement.com/elements/vector-lagrange.html}
   */
  template <>
  class P1Element<Math::Vector> final : public FiniteElementBase<P1Element<Math::Vector>>
  {
    using G = Geometry::Polytope::Type;

    public:
      /// Parent class
      using Parent = FiniteElementBase<P1Element>;

      /// Type of range
      using RangeType = Math::Vector;

      class LinearForm
      {
        public:
          LinearForm() = default;

          constexpr
          LinearForm(size_t vdim, size_t i, Geometry::Polytope::Type g)
            : m_vdim(vdim), m_i(i), m_g(g)
          {
            assert(m_vdim > 0);
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
            return v(s_nodes[m_vdim][m_g].col(m_i)).coeff(static_cast<size_t>(m_i % m_vdim));
          }

        private:
          size_t m_vdim;
          size_t m_i;
          Geometry::Polytope::Type m_g;
      };

      class BasisFunction
      {
        public:
          BasisFunction() = default;

          constexpr
          BasisFunction(size_t vdim, size_t i, Geometry::Polytope::Type g)
            : m_vdim(vdim), m_i(i), m_g(g)
          {
            assert(m_vdim > 0);
          }

          constexpr
          BasisFunction(const BasisFunction&) = default;

          constexpr
          BasisFunction(BasisFunction&&) = default;

          constexpr
          BasisFunction& operator=(const BasisFunction&) = default;

          constexpr
          BasisFunction& operator=(BasisFunction&&) = default;

          Math::Vector operator()(const Math::SpatialVector& r) const
          {
            Math::Vector res;
            operator()(res, r);
            return res;
          }

          void operator()(Math::Vector& out, const Math::SpatialVector& r) const;

        private:
          size_t m_vdim;
          size_t m_i;
          Geometry::Polytope::Type m_g;
      };

      class JacobianFunction
      {
        public:
          JacobianFunction() = default;

          constexpr
          JacobianFunction(size_t vdim, size_t i, Geometry::Polytope::Type g)
            : m_vdim(vdim), m_i(i), m_g(g)
          {
            assert(m_vdim > 0);
          }

          constexpr
          JacobianFunction(const JacobianFunction&) = default;

          constexpr
          JacobianFunction(JacobianFunction&&) = default;

          constexpr
          JacobianFunction& operator=(const JacobianFunction&) = default;

          constexpr
          JacobianFunction& operator=(JacobianFunction&&) = default;

          Math::Matrix operator()(const Math::SpatialVector& rc) const
          {
            Math::Matrix res;
            operator()(res, rc);
            return res;
          }

          void operator()(Math::Matrix& out, const Math::SpatialVector& rc) const;

        private:
          size_t m_vdim;
          size_t m_i;
          Geometry::Polytope::Type m_g;
      };

      P1Element() = default;

      P1Element(size_t vdim, Geometry::Polytope::Type geometry)
        : Parent(geometry),
          m_vdim(vdim)
      {}

      P1Element(const P1Element& other)
        : Parent(other),
          m_vdim(other.m_vdim)
      {}

      P1Element(P1Element&& other)
        : Parent(std::move(other)),
          m_vdim(other.m_vdim)
      {}

      P1Element& operator=(P1Element&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_vdim = std::move(other.m_vdim);
        }
        return *this;
      }

      P1Element& operator=(const P1Element& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_vdim = other.m_vdim;
        }
        return *this;
      }

      inline
      constexpr
      size_t getCount() const
      {
        return Geometry::Polytope::getVertexCount(getGeometry()) * m_vdim;
      }

      inline
      const Math::Matrix& getNodes() const
      {
        return s_nodes[m_vdim][getGeometry()];
      }

      inline
      const auto& getLinearForm(size_t i) const
      {
        assert(i < getCount());
        return s_ls[m_vdim][getGeometry()][i];
      }

      inline
      const auto& getBasis(size_t i) const
      {
        assert(i < getCount());
        return s_basis[m_vdim][getGeometry()][i];
      }

      inline
      const auto& getJacobian(size_t i) const
      {
        assert(i < getCount());
        return s_jacobian[m_vdim][getGeometry()][i];
      }

      inline
      size_t getOrder() const
      {
        switch (getGeometry())
        {
          case Geometry::Polytope::Type::Point:
            return 0;
          case Geometry::Polytope::Type::Segment:
          case Geometry::Polytope::Type::Triangle:
          case Geometry::Polytope::Type::Tetrahedron:
            return 1;
          case Geometry::Polytope::Type::Quadrilateral:
            return 2;
        }
        assert(false);
        return 0;
      }

    private:
      static const
      std::array<Geometry::GeometryIndexed<Math::Matrix>, RODIN_P1_MAX_VECTOR_DIMENSION> s_nodes;

      static const
      std::array<Geometry::GeometryIndexed<std::vector<LinearForm>>, RODIN_P1_MAX_VECTOR_DIMENSION> s_ls;

      static const
      std::array<Geometry::GeometryIndexed<std::vector<BasisFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION> s_basis;

      static const
      std::array<Geometry::GeometryIndexed<std::vector<JacobianFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION> s_jacobian;

      size_t m_vdim;
  };

  namespace Internal
  {
    std::array<Geometry::GeometryIndexed<Math::Matrix>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1Nodes();

    std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::LinearForm>>, RODIN_P1_MAX_VECTOR_DIMENSION> initVectorP1LinearForms();

    std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::BasisFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1Basis();

    std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::JacobianFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1Jacobian();
  }
}

#endif
