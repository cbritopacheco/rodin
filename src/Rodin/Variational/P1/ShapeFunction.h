#ifndef RODIN_VARIATIONAL_P1_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_P1_SHAPEFUNCTION_H

#include <utility>
#include <vector>

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/P1/ForwardDecls.h"
#include "Rodin/Variational/P1/P1.h"
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

      struct Cache
      {
        struct Key
        {
          Geometry::Polytope::Type geom = Geometry::Polytope::Type::Point;
          const QF::QuadratureFormulaBase* qf = nullptr;
          size_t qp = 0;
          bool valid = false;

          explicit operator bool() const noexcept { return valid; }

          bool operator==(const Key& other) const noexcept
          {
            if (!valid || !other.valid)
              return false;
            return geom == other.geom && qf == other.qf && qp == other.qp;
          }

          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            qf = nullptr;
            qp = 0;
            geom = Geometry::Polytope::Type::Point;
          }
        };

        std::vector<RangeType> basis;
        Key key;
      };

      ShapeFunction() = delete;

      constexpr
      ShapeFunction(const FESType& fes)
        : Parent(fes)
      {}

      constexpr
      ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_ip(nullptr),
          m_cache(other.m_cache)
      {}

      constexpr
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_cache(std::move(other.m_cache))
      {}

      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        if constexpr (std::is_same_v<RangeType, ScalarType>)
        {
          return P1Element<ScalarType>(polytope.getGeometry()).getCount();
        }
        else
        {
          static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);
          return P1Element<RangeType>(polytope.getGeometry(), this->getFiniteElementSpace().getVectorDimension()).getCount();
        }
      }

      constexpr
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /**
       * Fallback path (non-quadrature evaluations).
       * Keeps your previous behavior (pushforward of each basis at p).
       */
      ShapeFunction& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& p = ip.getPoint();
        const auto& qf = ip.getQuadratureFormula();
        const size_t qp = ip.getIndex();
        const auto& polytope = p.getPolytope();
        const Geometry::Polytope::Type geom = polytope.getGeometry();

        typename Cache::Key key;
        key.geom  = geom;
        key.qf    = &qf;
        key.qp    = qp;
        key.valid = true;

        const bool recompute = !(m_cache.key == key);

        if (recompute)
        {
          m_cache.key = key;

          if constexpr (std::is_same_v<RangeType, ScalarType>)
          {
            const P1Element<RangeType> fe(geom);
            const size_t ndof = fe.getCount();
            m_cache.basis.resize(ndof);
            const auto& rq = qf.getPoint(qp);
            for (size_t a = 0; a < ndof; ++a)
              m_cache.basis[a] = fe.getBasis(a)(rq);
          }
          else
          {
            static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);
            const P1Element<RangeType> fe(geom, this->getFiniteElementSpace().getVectorDimension());
            const size_t ndof = fe.getCount();
            m_cache.basis.resize(ndof);
            const auto& rq = qf.getPoint(qp);
            for (size_t a = 0; a < ndof; ++a)
              m_cache.basis[a] = fe.getBasis(a)(rq);
          }
        }

        return *this;
      }

      constexpr
      const RangeType& getBasis(size_t local) const
      {
        assert(m_cache.key);
        assert(local < m_cache.basis.size());
        return m_cache.basis[local];
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
      const IntegrationPoint* m_ip;
      Cache m_cache;
  };
}

#endif

