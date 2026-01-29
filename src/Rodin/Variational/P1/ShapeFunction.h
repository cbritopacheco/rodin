#ifndef RODIN_VARIATIONAL_P1_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_P1_SHAPEFUNCTION_H

#include <utility>
#include <vector>

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Variational/P1/ForwardDecls.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/IntegrationPoint.h"

namespace Rodin::Variational
{
  template <class Derived, class Range, class Mesh, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, P1<Range, Mesh>, Space>
    : public ShapeFunctionBase<ShapeFunction<Derived, P1<Range, Mesh>, Space>, P1<Range, Mesh>, Space>
  {
    public:
      using FESType = P1<Range, Mesh>;
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

      ShapeFunction& setIntegrationPoint(const IntegrationPoint& ip)
      {
        assert(ip.getPoint());
        return this->setPoint(*ip.getPoint());
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
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        return P1Element<ScalarType>(poly.getGeometry()).getOrder();
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

