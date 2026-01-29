#ifndef RODIN_VARIATIONAL_H1_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_H1_SHAPEFUNCTION_H

#include <utility>
#include <vector>

#include "Rodin/Variational/H1/ForwardDecls.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/IntegrationPoint.h"

namespace Rodin::Variational
{
  template <class Derived, size_t K, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, H1<K, Scalar, Mesh>, Space>
    : public ShapeFunctionBase<ShapeFunction<Derived, H1<K, Scalar, Mesh>, Space>, H1<K, Scalar, Mesh>, Space>
  {
    public:
      using FESType = H1<K, Scalar, Mesh>;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      using RangeType  = typename FormLanguage::Traits<FESType>::RangeType;

      using Parent =
        ShapeFunctionBase<
          ShapeFunction<Derived, FESType, SpaceType>,
          FESType,
          SpaceType>;

      ShapeFunction() = delete;

      constexpr
      ShapeFunction(const FESType& fes)
        : Parent(fes),
          m_p(nullptr)
      {}

      constexpr
      ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_basis(other.m_basis),
          m_p(nullptr)
      {}

      constexpr
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_basis(std::move(other.m_basis)),
          m_p(std::exchange(other.m_p, nullptr))
      {}

      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        const size_t d = polytope.getDimension();
        const Index  i = polytope.getIndex();
        return this->getFiniteElementSpace().getFiniteElement(d, i).getCount();
      }

      constexpr
      const Geometry::Point& getPoint() const
      {
        assert(m_p);
        return *m_p;
      }

      /**
       * Fallback path (non-quadrature evaluations).
       * Keeps your previous behavior (pushforward of each basis at p).
       */
      ShapeFunction& setPoint(const Geometry::Point& p)
      {
        if (m_p == &p)
          return *this;

        m_p = &p;

        const auto& polytope = p.getPolytope();
        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        const auto& fes = this->getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, idx);

        const size_t ndof = fe.getCount();
        m_basis.resize(ndof);

        for (size_t a = 0; a < ndof; ++a)
          m_basis[a] = fes.getPushforward({ d, idx }, fe.getBasis(a))(p);

        return *this;
      }

      /**
       * Fast path: quadrature evaluation via Tabulation.
       */
      ShapeFunction& setIntegrationPoint(const IntegrationPoint& ip)
      {
        const auto* pp = ip.getPoint();
        assert(pp);
        const Geometry::Point& p = *pp;

        // If you want: keep the same pointer short-circuit
        if (m_p == &p)
          return *this;

        const auto* qf = ip.getQuadratureFormula();
        assert(qf);

        const size_t qp = ip.getIndex();

        m_p = &p;

        const auto& polytope = p.getPolytope();
        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        const auto& fes = this->getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, idx);

        const size_t ndof = fe.getCount();
        m_basis.resize(ndof);

        const auto& tab = fe.getTabulation(*qf);

        for (size_t a = 0; a < ndof; ++a)
          m_basis[a] = tab.getBasis(qp, a);

        return *this;
      }

      constexpr
      const RangeType& getBasis(size_t local) const
      {
        assert(local < m_basis.size());
        return m_basis[local];
      }

      constexpr
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        return H1Element<K, Scalar>(geom.getGeometry()).getOrder();
      }

      virtual ShapeFunction* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::vector<RangeType> m_basis;
      const Geometry::Point* m_p;
  };

  template <class Derived, size_t K, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, H1<K, Math::Vector<Scalar>, Mesh>, Space>
    : public ShapeFunctionBase<
        ShapeFunction<Derived, H1<K, Math::Vector<Scalar>, Mesh>, Space>,
        H1<K, Math::Vector<Scalar>, Mesh>,
        Space>
  {
    public:
      using FESType = H1<K, Math::Vector<Scalar>, Mesh>;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType; // should be Scalar
      using RangeType  = typename FormLanguage::Traits<FESType>::RangeType;  // should be Math::Vector<Scalar>

      using Parent =
        ShapeFunctionBase<
          ShapeFunction<Derived, FESType, SpaceType>,
          FESType,
          SpaceType>;

      static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);

      ShapeFunction() = delete;

      constexpr
      ShapeFunction(const FESType& fes)
        : Parent(fes),
          m_p(nullptr)
      {}

      constexpr
      ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_basis(other.m_basis),
          m_p(nullptr)
      {}

      constexpr
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_basis(std::move(other.m_basis)),
          m_p(std::exchange(other.m_p, nullptr))
      {}

      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        const size_t d = polytope.getDimension();
        const Index  i = polytope.getIndex();
        return this->getFiniteElementSpace().getFiniteElement(d, i).getCount();
      }

      constexpr
      const Geometry::Point& getPoint() const
      {
        assert(m_p);
        return *m_p;
      }

      /**
       * Slow path: pushforward each basis at p.
       */
      ShapeFunction& setPoint(const Geometry::Point& p)
      {
        if (m_p == &p)
          return *this;

        m_p = &p;

        const auto& polytope = p.getPolytope();
        const size_t d   = polytope.getDimension();
        const Index  idx = polytope.getIndex();

        const auto& fes = this->getFiniteElementSpace();
        const auto& fe  = fes.getFiniteElement(d, idx);

        const size_t ndof = fe.getCount();
        m_basis.resize(ndof);

        for (size_t a = 0; a < ndof; ++a)
          m_basis[a] = fes.getPushforward({ d, idx }, fe.getBasis(a))(p);

        return *this;
      }

      /**
       * Intentionally slow: reuse setPoint (pushforward) even under quadrature.
       */
      ShapeFunction& setIntegrationPoint(const IntegrationPoint& ip)
      {
        const auto* pp = ip.getPoint();
        assert(pp);
        return setPoint(*pp);
      }

      constexpr
      const RangeType& getBasis(size_t local) const
      {
        assert(local < m_basis.size());
        return m_basis[local];
      }

      constexpr
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        return H1Element<K, Scalar>(geom.getGeometry()).getOrder();
      }

      virtual ShapeFunction* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::vector<RangeType> m_basis;
      const Geometry::Point* m_p;
  };
}

#endif
