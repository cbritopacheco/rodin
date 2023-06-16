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
    using G = Geometry::Polytope::Geometry;

    public:
      /// Type of range
      using RangeType = Scalar;

      /// Parent class
      using Parent = FiniteElementBase<P1Element<Scalar>>;

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

      class BasisFunction
      {
        public:
          constexpr
          BasisFunction(size_t i, Geometry::Polytope::Geometry g)
            : m_i(i), m_g(g)
          {
            assert(i < Geometry::Polytope::getVertexCount(g));
          }

          constexpr
          BasisFunction(const BasisFunction&) = default;

          constexpr
          BasisFunction& operator=(const BasisFunction&) = default;

          constexpr
          BasisFunction& operator=(BasisFunction&&) = default;

          inline
          constexpr
          Scalar operator()(const Math::Vector& r) const
          {
            switch (m_g)
            {
              case Geometry::Polytope::Geometry::Point:
                return 1;
              case Geometry::Polytope::Geometry::Segment:
              {
                switch (m_i)
                {
                  case 0:
                    return 1 - r.x();
                  case 1:
                    return r.x();
                  default:
                  {
                    assert(false);
                    return NAN;
                  }
                }
              }
              case Geometry::Polytope::Geometry::Triangle:
              {
                switch (m_i)
                {
                  case 0:
                    return -r.x() - r.y() + 1;
                  case 1:
                    return r.x();
                  case 2:
                    return r.y();
                  default:
                  {
                    assert(false);
                    return NAN;
                  }
                }
              }
              case Geometry::Polytope::Geometry::Quadrilateral:
              {
                switch (m_i)
                {
                  case 0:
                  {
                    auto x = r.x();
                    auto y = r.y();
                    return x * y - x - y + 1;
                  }
                  case 1:
                    return r.x() * (1 - r.y());
                  case 2:
                    return r.y() * (1 - r.x());
                  case 3:
                    return r.x() * r.y();
                  default:
                  {
                    assert(false);
                    return NAN;
                  }
                }
              }
              case Geometry::Polytope::Geometry::Tetrahedron:
              {
                switch (m_i)
                {
                  case 0:
                    return r.x() - r.y() - r.z() + 1;
                  case 1:
                    return r.x();
                  case 2:
                    return r.y();
                  case 3:
                    return r.z();
                  default:
                  {
                    assert(false);
                    return NAN;
                  }
                }
              }
            }
            assert(false);
            return NAN;
          }

        private:
          size_t m_i;
          Geometry::Polytope::Geometry m_g;
      };

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

          inline
          Math::Vector operator()(const Math::Vector& r) const
          {
            switch (m_g)
            {
              case Geometry::Polytope::Geometry::Point:
                return Math::Vector{{0, 0}};
              case Geometry::Polytope::Geometry::Segment:
              {
                switch (m_i)
                {
                  case 0:
                    return Math::Vector{{-1}};
                  case 1:
                    return Math::Vector{{1}};
                  default:
                  {
                    assert(false);
                    return Math::Vector{{NAN}};
                  }
                }
              }
              case Geometry::Polytope::Geometry::Triangle:
              {
                switch (m_i)
                {
                  case 0:
                    return Math::Vector{{-1, -1}};
                  case 1:
                    return Math::Vector{{1, 0}};
                  case 2:
                    return Math::Vector{{0, 1}};
                  default:
                  {
                    assert(false);
                    return Math::Vector{{NAN, NAN}};
                  }
                }
              }
              case Geometry::Polytope::Geometry::Quadrilateral:
              {
                switch (m_i)
                {
                  case 0:
                    return Math::Vector{{r.y() - 1, r.x() - 1}};
                  case 1:
                    return Math::Vector{{1 - r.y(), -r.y()}};
                  case 2:
                    return Math::Vector{{1 - r.x(), -r.x()}};
                  case 3:
                    return Math::Vector{{r.y(), r.x()}};
                  default:
                  {
                    assert(false);
                    return Math::Vector{{NAN, NAN}};
                  }
                }
              }
              case Geometry::Polytope::Geometry::Tetrahedron:
              {
                switch (m_i)
                {
                  case 0:
                    return Math::Vector{{-1, -1, -1}};
                  case 1:
                    return Math::Vector{{1, 0, 0}};
                  case 2:
                    return Math::Vector{{0, 1, 0}};
                  case 3:
                    return Math::Vector{{0, 0, 1}};
                  default:
                  {
                    assert(false);
                    return Math::Vector{{NAN, NAN, NAN}};
                  }
                }
              }
            }
            assert(false);
            return Math::Vector{{}};
          }

        private:
          size_t m_i;
          Geometry::Polytope::Geometry m_g;
      };

      constexpr
      P1Element(Geometry::Polytope::Geometry geometry)
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
      constexpr
      const auto& getLinearForm(size_t i) const
      {
        assert(i < getCount());
        return s_ls[getGeometry()][i];
      }

      inline
      constexpr
      const auto& getBasis(size_t i) const
      {
        assert(i < getCount());
        return s_basis[getGeometry()][i];
      }

      inline
      Math::Vector getBasis(const Math::Vector& r) const
      {
        const size_t n = getCount();
        Math::Vector res(n);
        for (size_t i = 0; i < n; i++)
          res.coeffRef(i) = getBasis(i)(r);
        return res;
      }

      inline
      constexpr
      const auto& getGradient(size_t i) const
      {
        assert(i < getCount());
        return s_gradient[getGeometry()][i];
      }

      inline
      Math::Matrix getGradient(const Math::Vector& r) const
      {
        const size_t n = getCount();
        const size_t d = Geometry::Polytope::getGeometryDimension(getGeometry());
        Math::Matrix res(d, n);
        for (size_t i = 0; i < n; i++)
          res.col(i) = getGradient(i)(r);
        return res;
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
    using G = Geometry::Polytope::Geometry;

    public:
      class LinearForm
      {
        public:
          LinearForm() = default;

          constexpr
          LinearForm(size_t vdim, size_t i, Geometry::Polytope::Geometry g)
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
            return v(s_nodes[m_vdim][m_g].col(m_i)).coeff(static_cast<size_t>(m_i / m_vdim));
          }

        private:
          size_t m_vdim;
          size_t m_i;
          Geometry::Polytope::Geometry m_g;
      };

      class BasisFunction
      {
        public:
          BasisFunction() = default;

          constexpr
          BasisFunction(size_t vdim, size_t i, Geometry::Polytope::Geometry g)
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

          inline
          Math::Vector operator()(const Math::Vector& r) const
          {
            Math::Vector res = Math::Vector::Zero(m_vdim);
            const size_t i = m_i % m_vdim;
            const size_t k = m_i / m_vdim;
            assert(k < Geometry::Polytope::getVertexCount(m_g));
            switch (m_g)
            {
              case Geometry::Polytope::Geometry::Point:
              {
                res.coeffRef(0) = 1;
                return res;
              }
              case Geometry::Polytope::Geometry::Segment:
              {
                switch (k)
                {
                  case 0:
                  {
                    res.coeffRef(i) = 1 - r.x();
                    return res;
                  }
                  case 1:
                  {
                    res.coeffRef(i) = r.x();
                    return res;
                  }
                  default:
                  {
                    assert(false);
                    res.setConstant(NAN);
                    break;
                  }
                }
              }
              case Geometry::Polytope::Geometry::Triangle:
              {
                switch (k)
                {
                  case 0:
                  {
                    res(i) = -r.x() - r.y() + 1;
                    return res;
                  }
                  case 1:
                  {
                    res(i) = r.x();
                    return res;
                  }
                  case 2:
                  {
                    res(i) = r.y();
                    return res;
                  }
                  default:
                  {
                    assert(false);
                    res.setConstant(NAN);
                    break;
                  }
                }
              }
              case Geometry::Polytope::Geometry::Quadrilateral:
              {
                switch (k)
                {
                  case 0:
                  {
                    res(i) = r.x() * r.y() - r.x() - r.y() + 1;
                    return res;
                  }
                  case 1:
                  {
                    res(i) = r.x() * (1 - r.y());
                    return res;
                  }
                  case 2:
                  {
                    res(i) = r.y() * (1 - r.x());
                    return res;
                  }
                  case 3:
                  {
                    res(i) = r.x() * r.y();
                    return res;
                  }
                  default:
                  {
                    assert(false);
                    res.setConstant(NAN);
                    break;
                  }
                }
              }
              case Geometry::Polytope::Geometry::Tetrahedron:
              {
                switch (k)
                {
                  case 0:
                  {
                    res(i) = -r.x() - r.y() - r.z() + 1;
                    return res;
                  }
                  case 1:
                  {
                    res(i) = r.x();
                    return res;
                  }
                  case 2:
                  {
                    res(i) = r.y();
                    return res;
                  }
                  case 3:
                  {
                    res(i) = r.z();
                    return res;
                  }
                  default:
                  {
                    assert(false);
                    res.setConstant(NAN);
                    break;
                  }
                }
              }
            }
            assert(false);
            res.setConstant(NAN);
            return res;
          }

        private:
          size_t m_vdim;
          size_t m_i;
          Geometry::Polytope::Geometry m_g;
      };

      class JacobianFunction
      {
        public:
          JacobianFunction() = default;

          constexpr
          JacobianFunction(size_t vdim, size_t i, Geometry::Polytope::Geometry g)
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

          inline
          Math::Matrix operator()(const Math::Vector& rc) const
          {
            assert(false);
            switch (m_g)
            {
              default:
                return Math::Matrix{{m_i * 1.0}};
            }
            return {};
          }

        private:
          size_t m_vdim;
          size_t m_i;
          Geometry::Polytope::Geometry m_g;
      };

      /// Type of range
      using RangeType = Math::Vector;

      /// Parent class
      using Parent = FiniteElementBase<P1Element>;

      constexpr
      P1Element(size_t vdim, Geometry::Polytope::Geometry geometry)
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

    private:
      static const
      std::array<Geometry::GeometryIndexed<Math::Matrix>, RODIN_P1_MAX_VECTOR_DIMENSION> s_nodes;

      static const
      std::array<Geometry::GeometryIndexed<std::vector<LinearForm>>, RODIN_P1_MAX_VECTOR_DIMENSION> s_ls;

      static const
      std::array<Geometry::GeometryIndexed<std::vector<BasisFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION> s_basis;

      static const
      std::array<Geometry::GeometryIndexed<std::vector<JacobianFunction>>, RODIN_P1_MAX_VECTOR_DIMENSION> s_jacobian;

      const size_t m_vdim;
  };

  /**
   * @brief Alias for P1Element<Scalar>
   */
  using ScalarP1Element = P1Element<Scalar>;

  /**
   * @brief Alias for P1Element<Math::Vector>
   */
  using VectorP1Element = P1Element<Math::Vector>;

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
