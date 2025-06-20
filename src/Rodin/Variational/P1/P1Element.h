/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_P1ELEMENT_H
#define RODIN_VARIATIONAL_P1_P1ELEMENT_H

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
  class P1Element<Real> final : public FiniteElementBase<P1Element<Real>>
  {
    using G = Geometry::Polytope::Type;

    public:
      /// Parent class
      using Parent = FiniteElementBase<P1Element<Real>>;

      /// Type of range
      using RangeType = Real;

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

          template <class T>
          constexpr
          auto operator()(const T& v) const
          {
            return v(s_nodes[m_g].col(m_i));
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
          using ReturnType = Real;

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
          {
            assert(local < Geometry::Polytope::getVertexCount(g));
          }

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

          constexpr
          DerivativeFunction<1> getDerivative(size_t i) const
          {
            return DerivativeFunction<1>(i, m_local, m_g);
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
        return Geometry::Polytope::getVertexCount(getGeometry());
      }

      constexpr
      const Math::PointMatrix& getNodes() const
      {
        return s_nodes[getGeometry()];
      }

      constexpr
      auto getLinearForm(size_t i) const
      {
        return LinearForm(i, getGeometry());
      }

      constexpr
      auto getBasis(size_t i) const
      {
        return BasisFunction(i, getGeometry());
      }

      constexpr
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
          case Geometry::Polytope::Type::Wedge:
            return 2;
        }
        assert(false);
        return 0;
      }

    private:
      static const Geometry::GeometryIndexed<Math::PointMatrix> s_nodes;
  };

  template <>
  class P1Element<Complex> final : public FiniteElementBase<P1Element<Complex>>
  {
    using G = Geometry::Polytope::Type;

    public:
      /// Parent class
      using Parent = FiniteElementBase<P1Element<Complex>>;

      /// Type of range
      using RangeType = Complex;

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
          using ReturnType = Complex;

          template <size_t Order>
          class DerivativeFunction
          {
            public:
              using ReturnType = Complex;

              constexpr
              DerivativeFunction(size_t i, size_t local, Geometry::Polytope::Type g)
                : m_derivative(i, local, g)
              {}

              constexpr
              DerivativeFunction(const DerivativeFunction&) = default;

              constexpr
              DerivativeFunction& operator=(const DerivativeFunction&) = default;

              constexpr
              DerivativeFunction& operator=(DerivativeFunction&&) = default;

              constexpr
              void operator()(ReturnType& out, const Math::SpatialVector<Real>& r) const
              {
                out.real(m_derivative(r));
                out.imag(0);
              }

              constexpr
              ReturnType operator()(const Math::SpatialVector<Real>& r) const
              {
                return Complex(m_derivative(r), 0);
              }

            private:
              RealP1Element::BasisFunction::DerivativeFunction<Order> m_derivative;
          };

          constexpr
          BasisFunction(size_t local, Geometry::Polytope::Type g)
            : m_local(local), m_g(g),
              m_basis(local, g)
          {}

          constexpr
          BasisFunction(const BasisFunction&) = default;

          constexpr
          ReturnType operator()(const Math::SpatialVector<Real>& r) const
          {
            return Complex(m_basis(r), 0);
          }

          constexpr
          void operator()(ReturnType& out, const Math::SpatialVector<Real>& r) const
          {
            out.real(m_basis(r));
            out.imag(0);
          }

          template <size_t Order>
          constexpr
          DerivativeFunction<Order> getDerivative(size_t i) const
          {
            return DerivativeFunction<Order>(i, m_local, m_g);
          }

          constexpr
          DerivativeFunction<1> getDerivative(size_t i) const
          {
            return DerivativeFunction<1>(i, m_local, m_g);
          }

        private:
          const size_t m_local;
          const Geometry::Polytope::Type m_g;
          const RealP1Element::BasisFunction m_basis;
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
      constexpr
      size_t getCount() const
      {
        return Geometry::Polytope::getVertexCount(getGeometry());
      }

      constexpr
      const Math::PointMatrix& getNodes() const
      {
        return s_nodes[getGeometry()];
      }

      const auto& getLinearForm(size_t local) const
      {
        assert(local < getCount());
        return s_ls[getGeometry()][local];
      }

      const BasisFunction& getBasis(size_t local) const
      {
        assert(local < getCount());
        return s_basis[getGeometry()][local];
      }

      constexpr
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
          case Geometry::Polytope::Type::Wedge:
            return 2;
        }
        assert(false);
        return 0;
      }

    private:
      static const Geometry::GeometryIndexed<Math::PointMatrix> s_nodes;
      static const Geometry::GeometryIndexed<std::vector<LinearForm>> s_ls;
      static const Geometry::GeometryIndexed<std::vector<BasisFunction>> s_basis;
  };

  /**
   * @ingroup FiniteElements
   * @ingroup P1ElementSpecializations
   * @brief Degree 1 vector Lagrange element
   *
   * @see @m_defelement{Vector Lagrange,https://defelement.com/elements/vector-lagrange.html}
   */
  template <>
  class P1Element<Math::Vector<Real>> final : public FiniteElementBase<P1Element<Math::Vector<Real>>>
  {
    using G = Geometry::Polytope::Type;

    public:
      /// Parent class
      using Parent = FiniteElementBase<P1Element>;

      /// Type of range
      using RangeType = Math::Vector<Real>;

      class LinearForm
      {
        public:
          constexpr
          LinearForm()
            : m_vdim(0), m_i(0), m_g(Geometry::Polytope::Type::Point)
          {}

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
          constexpr
          BasisFunction()
            : m_vdim(0), m_i(0), m_g(Geometry::Polytope::Type::Point)
          {}

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

          Math::Vector<Real> operator()(const Math::SpatialVector<Real>& r) const
          {
            Math::Vector<Real> res;
            operator()(res, r);
            return res;
          }

          void operator()(Math::Vector<Real>& out, const Math::SpatialVector<Real>& r) const;

        private:
          size_t m_vdim;
          size_t m_i;
          Geometry::Polytope::Type m_g;
      };

      P1Element() = default;

      constexpr
      P1Element(size_t vdim, Geometry::Polytope::Type geometry)
        : Parent(geometry),
          m_vdim(vdim)
      {}

      constexpr
      P1Element(const P1Element& other)
        : Parent(other),
          m_vdim(other.m_vdim)
      {}

      constexpr
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

      constexpr
      size_t getCount() const
      {
        return Geometry::Polytope::getVertexCount(getGeometry()) * m_vdim;
      }

      constexpr
      const Math::PointMatrix& getNodes() const
      {
        return s_nodes[m_vdim][getGeometry()];
      }

      constexpr
      auto getLinearForm(size_t i) const
      {
        assert(i < getCount());
        return s_ls[m_vdim][getGeometry()][i];
      }

      constexpr
      auto getBasis(size_t i) const
      {
        assert(i < getCount());

      }


      constexpr
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
          case Geometry::Polytope::Type::Wedge:
            return 2;
        }
        assert(false);
        return 0;
      }

    private:
      static const
      std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_P1_MAX_VECTOR_DIMENSION> s_nodes;

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
    std::array<Geometry::GeometryIndexed<Math::PointMatrix>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1Nodes();

    std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::LinearForm>>, RODIN_P1_MAX_VECTOR_DIMENSION> initVectorP1LinearForms();

    std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::BasisFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1Basis();

    std::array<Geometry::GeometryIndexed<std::vector<VectorP1Element::JacobianFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION>
    initVectorP1Jacobian();
  }
}

#include "P1Element.hpp"

#endif
