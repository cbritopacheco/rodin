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
   * @brief Degree 0 scalar Lagrange element
   * @see @m_defelement{Lagrange,https://defelement.com/elements/lagrange.html}
   */
  template <>
  class P0Element<Scalar> final : public FiniteElementBase<P0Element<Scalar>>
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
          LinearForm()
            : LinearForm(Geometry::Polytope::Type::Point)
          {}

          constexpr
          LinearForm(Geometry::Polytope::Type g)
            : m_g(g)
          {}

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
            return v(s_nodes[m_g].col(0));
          }

        private:
          Geometry::Polytope::Type m_g;
      };

      /**
       * @brief Represents a basis function of a P0 scalar element.
       */
      class BasisFunction
      {
        public:
          using ReturnType = Scalar;

          constexpr
          BasisFunction()
          {}

          constexpr
          BasisFunction(Geometry::Polytope::Type)
          {}

          constexpr
          BasisFunction(const BasisFunction&) = default;

          constexpr
          BasisFunction& operator=(const BasisFunction&) = default;

          constexpr
          BasisFunction& operator=(BasisFunction&&) = default;

          inline
          constexpr
          Scalar operator()(const Math::SpatialVector& r) const
          {
            return 1;
          }
      };

      /**
       * @brief Represents a gradient basis function of a P0 scalar element.
       */
      class GradientFunction
      {
        public:
          constexpr
          GradientFunction()
            : GradientFunction(Geometry::Polytope::Type::Point)
          {}

          constexpr
          GradientFunction(Geometry::Polytope::Type g)
            : m_g(g)
          {}

          constexpr
          GradientFunction(const GradientFunction&) = default;

          constexpr
          GradientFunction& operator=(const GradientFunction&) = default;

          constexpr
          GradientFunction& operator=(GradientFunction&&) = default;

          inline
          auto operator()(const Math::SpatialVector& r) const
          {
            return Math::SpatialVector::Zero(Geometry::Polytope::getGeometryDimension(m_g));
          }

          inline
          constexpr
          void operator()(Math::SpatialVector& out, const Math::SpatialVector& r) const
          {
            assert(r.size() == Geometry::Polytope::getGeometryDimension(m_g));
            out.resize(Geometry::Polytope::getGeometryDimension(m_g));
            out.setZero();
          }

        private:
          Geometry::Polytope::Type m_g;
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
      constexpr
      size_t getCount() const
      {
        return 1;
      }

      inline
      constexpr
      const Math::PointMatrix& getNodes() const
      {
        return s_nodes[getGeometry()];
      }

      inline
      constexpr
      const auto& getLinearForm(size_t i) const
      {
        assert(i == 0);
        return s_ls[getGeometry()];
      }

      inline
      constexpr
      const auto& getBasis(size_t i) const
      {
        assert(i == 0);
        return s_basis[getGeometry()];
      }

      inline
      constexpr
      const auto& getGradient(size_t i) const
      {
        assert(i == 0);
        return s_gradient[getGeometry()];
      }

    private:
      static const Geometry::GeometryIndexed<Math::PointMatrix> s_nodes;
      static const Geometry::GeometryIndexed<LinearForm> s_ls;
      static const Geometry::GeometryIndexed<BasisFunction> s_basis;
      static const Geometry::GeometryIndexed<GradientFunction> s_gradient;
  };

  // /**
  //  * @ingroup FiniteElements
  //  * @ingroup P0ElementSpecializations
  //  * @brief Degree 0 vector Lagrange element
  //  *
  //  * @see @m_defelement{Vector Lagrange,https://defelement.com/elements/vector-lagrange.html}
  //  */
  // template <>
  // class P0Element<Math::Vector> final : public FiniteElementBase<P0Element<Math::Vector>>
  // {
  //   using G = Geometry::Polytope::Type;

  //   public:
  //     /// Parent class
  //     using Parent = FiniteElementBase<P0Element>;

  //     /// Type of range
  //     using RangeType = Math::Vector;

  //     class LinearForm
  //     {
  //       public:
  //         constexpr
  //         LinearForm()
  //           : m_vdim(0), m_i(0), m_g(Geometry::Polytope::Type::Point)
  //         {}

  //         constexpr
  //         LinearForm(size_t vdim, size_t i, Geometry::Polytope::Type g)
  //           : m_vdim(vdim), m_i(i), m_g(g)
  //         {
  //           assert(m_vdim > 0);
  //         }

  //         constexpr
  //         LinearForm(const LinearForm&) = default;

  //         constexpr
  //         LinearForm(LinearForm&&) = default;

  //         constexpr
  //         LinearForm& operator=(const LinearForm&) = default;

  //         constexpr
  //         LinearForm& operator=(LinearForm&&) = default;

  //         template <class T>
  //         inline
  //         constexpr
  //         auto operator()(const T& v) const
  //         {
  //           return v(s_nodes[m_vdim][m_g].col(m_i)).coeff(static_cast<size_t>(m_i % m_vdim));
  //         }

  //       private:
  //         size_t m_vdim;
  //         size_t m_i;
  //         Geometry::Polytope::Type m_g;
  //     };

  //     class BasisFunction
  //     {
  //       public:
  //         constexpr
  //         BasisFunction()
  //           : m_vdim(0), m_i(0), m_g(Geometry::Polytope::Type::Point)
  //         {}

  //         constexpr
  //         BasisFunction(size_t vdim, size_t i, Geometry::Polytope::Type g)
  //           : m_vdim(vdim), m_i(i), m_g(g)
  //         {
  //           assert(m_vdim > 0);
  //         }

  //         constexpr
  //         BasisFunction(const BasisFunction&) = default;

  //         constexpr
  //         BasisFunction(BasisFunction&&) = default;

  //         constexpr
  //         BasisFunction& operator=(const BasisFunction&) = default;

  //         constexpr
  //         BasisFunction& operator=(BasisFunction&&) = default;

  //         inline
  //         Math::Vector operator()(const Math::SpatialVector& r) const
  //         {
  //           Math::Vector res;
  //           operator()(res, r);
  //           return res;
  //         }

  //         inline
  //         void operator()(Math::Vector& out, const Math::SpatialVector& r) const
  //         {
  //         }

  //       private:
  //         size_t m_vdim;
  //         size_t m_i;
  //         Geometry::Polytope::Type m_g;
  //     };

  //     class JacobianFunction
  //     {
  //       public:
  //         constexpr
  //         JacobianFunction()
  //           : m_vdim(0), m_g(Geometry::Polytope::Type::Point)
  //         {}

  //         constexpr
  //         JacobianFunction(size_t vdim, size_t i, Geometry::Polytope::Type g)
  //           : m_vdim(vdim), m_g(g)
  //         {
  //           assert(m_vdim > 0);
  //         }

  //         constexpr
  //         JacobianFunction(const JacobianFunction&) = default;

  //         constexpr
  //         JacobianFunction(JacobianFunction&&) = default;

  //         constexpr
  //         JacobianFunction& operator=(const JacobianFunction&) = default;

  //         constexpr
  //         JacobianFunction& operator=(JacobianFunction&&) = default;

  //         inline
  //         auto operator()(const Math::SpatialVector& rc) const
  //         {
  //           return Math::SpatialMatrix::Zero(m_vdim, Geometry::Polytope::getGeometryDimension(m_g));
  //         }

  //         inline
  //         constexpr
  //         void operator()(Math::SpatialMatrix& out, const Math::SpatialVector& rc) const
  //         {
  //           out.resize(m_vdim, Geometry::Polytope::getGeometryDimension(m_g));
  //           out.setZero();
  //         }

  //       private:
  //         size_t m_vdim;
  //         Geometry::Polytope::Type m_g;
  //     };

  //     P0Element() = default;

  //     P0Element(size_t vdim, Geometry::Polytope::Type geometry)
  //       : Parent(geometry),
  //         m_vdim(vdim)
  //     {}

  //     P0Element(const P0Element& other)
  //       : Parent(other),
  //         m_vdim(other.m_vdim)
  //     {}

  //     P0Element(P0Element&& other)
  //       : Parent(std::move(other)),
  //         m_vdim(other.m_vdim)
  //     {}

  //     P0Element& operator=(P0Element&& other)
  //     {
  //       if (this != &other)
  //       {
  //         Parent::operator=(std::move(other));
  //         m_vdim = std::move(other.m_vdim);
  //       }
  //       return *this;
  //     }

  //     P0Element& operator=(const P0Element& other)
  //     {
  //       if (this != &other)
  //       {
  //         Parent::operator=(other);
  //         m_vdim = other.m_vdim;
  //       }
  //       return *this;
  //     }

  //     inline
  //     constexpr
  //     size_t getCount() const
  //     {
  //       return m_vdim;
  //     }

  //     inline
  //     constexpr
  //     const Math::PointMatrix& getNodes() const
  //     {
  //       return s_nodes[m_vdim][getGeometry()];
  //     }

  //     inline
  //     constexpr
  //     const auto& getLinearForm(size_t i) const
  //     {
  //       assert(i < getCount());
  //       return s_ls[m_vdim][getGeometry()];
  //     }

  //     inline
  //     constexpr
  //     const auto& getBasis(size_t i) const
  //     {
  //       assert(i < getCount());
  //       return s_basis[m_vdim][getGeometry()];
  //     }

  //     inline
  //     constexpr
  //     const auto& getJacobian(size_t i) const
  //     {
  //       assert(i < getCount());
  //       return s_jacobian[m_vdim][getGeometry()];
  //     }

  //   private:
  //     static const
  //     std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_P0_MAX_VECTOR_DIMENSION> s_nodes;

  //     static const
  //     std::array<Geometry::GeometryIndexed<LinearForm>, RODIN_P0_MAX_VECTOR_DIMENSION> s_ls;

  //     static const
  //     std::array<Geometry::GeometryIndexed<BasisFunction>, RODIN_P0_MAX_VECTOR_DIMENSION> s_basis;

  //     static const
  //     std::array<Geometry::GeometryIndexed<JacobianFunction>, RODIN_P0_MAX_VECTOR_DIMENSION> s_jacobian;

  //     size_t m_vdim;
  // };

  // namespace Internal
  // {
  //   std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_P0_MAX_VECTOR_DIMENSION>
  //   initVectorP0Nodes();
  // }
}

#endif
