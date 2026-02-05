/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0g_P0GELEMENT_H
#define RODIN_VARIATIONAL_P0g_P0GELEMENT_H

/**
 * @file
 * @brief P0g (global constant) finite element implementation.
 *
 * IMPORTANT:
 * - P0g is a *space* with 1 global DOF (scalar) or vdim global DOFs (vector),
 *   but its *element* is still the constant basis function on each cell.
 *
 * This element is therefore essentially identical to P0Element:
 * - Scalar basis:  phi(x) = 1
 * - Vector basis:  phi_i(x) = e_i (constant unit vectors)
 * - All derivatives are zero.
 *
 * Having a dedicated P0gElement is mostly for naming/clarity and to decouple
 * P0g from P0Element headers if you prefer.
 */

#include <cassert>
#include <utility>
#include <vector>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Variational/FiniteElement.h"

#include "ForwardDecls.h"

namespace Rodin::FormLanguage
{
  template <class Range>
  struct Traits<Variational::P0gElement<Range>>
  {
    using ScalarType = typename FormLanguage::Traits<Range>::ScalarType;
    using RangeType = Range;
  };
}

namespace Rodin::Variational
{
  // ----------------------------------------------------------------------------
  // Scalar P0gElement
  // ----------------------------------------------------------------------------
  template <class Scalar>
  class P0gElement final : public FiniteElementBase<P0gElement<Scalar>>
  {
    using G = Geometry::Polytope::Type;

  public:
    using Parent = FiniteElementBase<P0gElement<Scalar>>;
    using ScalarType = Scalar;
    using RangeType = Scalar;

    class LinearForm
    {
    public:
      constexpr explicit LinearForm(G g) : m_g(g) {}
      constexpr LinearForm(const LinearForm&) = default;

      template <class T>
      constexpr ScalarType operator()(const T& v) const
      {
        return v(P0gElement<ScalarType>(m_g).getNode(0));
      }

    private:
      const G m_g;
    };

    class BasisFunction
    {
    public:
      using ReturnType = Scalar;

      template <size_t Order>
      class DerivativeFunction
      {
      public:
        constexpr DerivativeFunction() = default;

        constexpr DerivativeFunction(const DerivativeFunction&) = default;

        constexpr ReturnType operator()(const Math::SpatialVector<Real>&) const
        {
          return ReturnType(0);
        }
      };

      constexpr BasisFunction() = default;

      constexpr BasisFunction(const BasisFunction&) = default;

      constexpr ReturnType operator()(const Math::SpatialVector<Real>&) const
      {
        return ReturnType(1);
      }

      template <size_t Order>
      constexpr DerivativeFunction<Order> getDerivative(size_t) const
      {
        return DerivativeFunction<Order>();
      }
    };

    constexpr P0gElement() = default;

    constexpr explicit P0gElement(G geometry)
      : Parent(geometry)
    {}

    constexpr P0gElement(const P0gElement&) = default;

    constexpr P0gElement(P0gElement&& other)
      : Parent(std::move(other))
    {}

    constexpr P0gElement& operator=(const P0gElement& other)
    {
      Parent::operator=(other);
      return *this;
    }

    constexpr ~P0gElement() override = default;

    constexpr size_t getCount() const { return 1; }

    const Math::SpatialVector<Real>& getNode(size_t i) const
    {
      assert(i == 0);
      switch (this->getGeometry())
      {
        case G::Point:
        {
          static thread_local const Math::SpatialVector<Real> s_node{};
          return s_node;
        }
        case G::Segment:
        {
          static thread_local const Math::SpatialVector<Real> s_node{ 0.5 };
          return s_node;
        }
        case G::Triangle:
        {
          static thread_local const Math::SpatialVector<Real> s_node{
            { Real(1) / Real(3), Real(1) / Real(3) }
          };
          return s_node;
        }
        case G::Quadrilateral:
        {
          static thread_local const Math::SpatialVector<Real> s_node{ { 0.5, 0.5 } };
          return s_node;
        }
        case G::Tetrahedron:
        {
          static thread_local const Math::SpatialVector<Real> s_node{ { 0.25, 0.25, 0.25 } };
          return s_node;
        }
        case G::Wedge:
        {
          static thread_local const Math::SpatialVector<Real> s_node{
            { Real(1) / Real(3), Real(1) / Real(3), 0.5 }
          };
          return s_node;
        }
        case G::Hexahedron:
        {
          static thread_local const Math::SpatialVector<Real> s_node{ { 0.5, 0.5, 0.5 } };
          return s_node;
        }
      }
      assert(false);
      static thread_local const Math::SpatialVector<Real> s_null{};
      return s_null;
    }

    constexpr LinearForm getLinearForm(size_t) const
    {
      return LinearForm(this->getGeometry());
    }

    constexpr BasisFunction getBasis(size_t) const
    {
      return BasisFunction();
    }

    constexpr size_t getOrder() const { return 0; }
  };

  // ----------------------------------------------------------------------------
  // Vector P0gElement
  // ----------------------------------------------------------------------------
  template <class Scalar>
  class P0gElement<Math::Vector<Scalar>> final
    : public FiniteElementBase<P0gElement<Math::Vector<Scalar>>>
  {
    using G = Geometry::Polytope::Type;

  public:
    using Parent = FiniteElementBase<P0gElement<Math::Vector<Scalar>>>;
    using ScalarType = Scalar;
    using RangeType = Math::Vector<Scalar>;

    class LinearForm
    {
    public:
      constexpr LinearForm()
        : m_vdim(0), m_local(0), m_g(G::Point)
      {}

      constexpr LinearForm(size_t vdim, size_t local, G g)
        : m_vdim(vdim), m_local(local), m_g(g)
      {}

      constexpr LinearForm(const LinearForm&) = default;

      template <class T>
      ScalarType operator()(const T& v) const
      {
        static thread_local RangeType s_out;
        const auto& xi = P0gElement<ScalarType>(m_g).getNode(m_local / m_vdim);
        s_out = v(xi);
        return s_out.coeff(m_local % m_vdim);
      }

    private:
      const size_t m_vdim;
      const size_t m_local;
      const G m_g;
    };

    class BasisFunction
    {
    public:
      using ReturnType = RangeType;

      template <size_t Order>
      class DerivativeFunction
      {
      public:
        constexpr DerivativeFunction(size_t, size_t, size_t, size_t, G) {}
        constexpr DerivativeFunction(const DerivativeFunction&) = default;

        constexpr ScalarType operator()(const Math::SpatialVector<Real>&) const
        {
          return ScalarType(0);
        }
      };

      constexpr BasisFunction()
        : m_vdim(0), m_local(0), m_g(G::Point)
      {}

      constexpr BasisFunction(size_t vdim, size_t local, G g)
        : m_vdim(vdim), m_local(local), m_g(g)
      {}

      constexpr BasisFunction(const BasisFunction&) = default;

      const ReturnType& operator()(const Math::SpatialVector<Real>&) const
      {
        static thread_local ReturnType s_out;
        s_out.resize(m_vdim);
        s_out.setZero();
        s_out.coeffRef(m_local % m_vdim) = ScalarType(1);
        return s_out;
      }

      template <size_t Order>
      constexpr DerivativeFunction<Order> getDerivative(size_t i, size_t j) const
      {
        return DerivativeFunction<Order>(i, j, m_vdim, m_local, m_g);
      }

    private:
      const size_t m_vdim;
      const size_t m_local;
      const G m_g;
    };

    P0gElement()
      : Parent(G::Point), m_vdim(0)
    {}

    /// Backward-compatible: vdim defaults to spatial dimension of geometry
    constexpr explicit P0gElement(G geometry)
      : P0gElement(geometry, Geometry::Polytope::Traits(geometry).getDimension())
    {}

    constexpr P0gElement(G geometry, size_t vdim)
      : Parent(geometry), m_vdim(vdim)
    {
      const size_t count = getCount();
      m_lfs.reserve(count);
      m_bs.reserve(count);
      for (size_t i = 0; i < count; ++i)
      {
        m_lfs.emplace_back(vdim, i, geometry);
        m_bs.emplace_back(vdim, i, geometry);
      }
    }

    constexpr P0gElement(const P0gElement& other)
      : Parent(other)
      , m_vdim(other.m_vdim)
      , m_lfs(other.m_lfs)
      , m_bs(other.m_bs)
    {}

    constexpr P0gElement(P0gElement&& other)
      : Parent(std::move(other))
      , m_vdim(std::exchange(other.m_vdim, 0))
      , m_lfs(std::move(other.m_lfs))
      , m_bs(std::move(other.m_bs))
    {}

    constexpr ~P0gElement() override = default;

    constexpr P0gElement& operator=(const P0gElement& other)
    {
      Parent::operator=(other);
      m_vdim = other.m_vdim;
      m_lfs  = other.m_lfs;
      m_bs   = other.m_bs;
      return *this;
    }

    constexpr P0gElement& operator=(P0gElement&& other)
    {
      Parent::operator=(std::move(other));
      m_vdim = std::exchange(other.m_vdim, 0);
      m_lfs  = std::move(other.m_lfs);
      m_bs   = std::move(other.m_bs);
      return *this;
    }

    constexpr size_t getCount() const
    {
      return m_vdim;
    }

    constexpr auto getLinearForm(size_t local) const
    {
      return m_lfs.at(local);
    }

    constexpr BasisFunction getBasis(size_t local) const
    {
      return m_bs.at(local);
    }

    constexpr const Math::SpatialVector<Real>& getNode(size_t local) const
    {
      // All components share the same barycentric node
      return P0gElement<ScalarType>(this->getGeometry()).getNode(local / m_vdim);
    }

    constexpr size_t getOrder() const { return 0; }

  private:
    size_t m_vdim;
    std::vector<LinearForm>    m_lfs;
    std::vector<BasisFunction> m_bs;
  };
}

#endif
