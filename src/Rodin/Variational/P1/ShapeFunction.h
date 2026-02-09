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
        struct StructureKey
        {
          Geometry::Polytope::Type geom = Geometry::Polytope::Type::Point;
          size_t vdim = 1;
          bool valid = false;

          explicit operator bool() const noexcept { return valid; }

          bool operator==(const StructureKey& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return geom == o.geom && vdim == o.vdim;
          }

          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            geom = Geometry::Polytope::Type::Point;
            vdim = 1;
          }
        };

        struct ValueKey
        {
          const QF::QuadratureFormulaBase* qf = nullptr;
          size_t qp = 0;
          bool valid = false;

          explicit operator bool() const noexcept { return valid; }

          bool operator==(const ValueKey& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return qf == o.qf && qp == o.qp;
          }

          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            qf = nullptr;
            qp = 0;
          }
        };

        // For scalar: size = nv
        // For vector: size = nv * vdim, where local = a*vdim + c
        std::vector<RangeType>  basis;

        // Scalar vertex basis values: size = nv
        std::vector<ScalarType> phi_vertex;

        StructureKey skey;
        ValueKey vkey;
      };

      ShapeFunction() = delete;

      constexpr
      ShapeFunction(const FESType& fes)
        : Parent(fes),
          m_ip(nullptr)
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

        const auto& p  = ip.getPoint();
        const auto& qf = ip.getQuadratureFormula();
        const size_t qp = ip.getIndex();

        const auto& poly = p.getPolytope();
        const auto geom  = poly.getGeometry();

        const size_t vdim =
          [] (const auto& self) -> size_t
          {
            if constexpr (std::is_same_v<RangeType, ScalarType>)
              return 1;
            else
              return self.getFiniteElementSpace().getVectorDimension();
          }(*this);

        // ---- structure cache: allocate/size once per (geom, vdim)
        typename Cache::StructureKey skey;
        skey.geom  = geom;
        skey.vdim  = vdim;
        skey.valid = true;

        const bool structure_changed = !(m_cache.skey == skey);
        if (structure_changed)
        {
          m_cache.skey = skey;
          m_cache.vkey = {}; // invalidate value cache

          const size_t nv = Geometry::Polytope::Traits(geom).getVertexCount();
          m_cache.phi_vertex.resize(nv);

          if constexpr (std::is_same_v<RangeType, ScalarType>)
          {
            m_cache.basis.resize(nv);
          }
          else
          {
            static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);
            const size_t ndof = nv * vdim;
            m_cache.basis.resize(ndof);

            // Initialize each stored vector once; later we only touch one component.
            for (auto& b : m_cache.basis)
            {
              b.resize(vdim);
              b.setZero();
            }
          }
        }

        // ---- value cache: update once per (qf, qp)
        typename Cache::ValueKey vkey;
        vkey.qf    = &qf;
        vkey.qp    = qp;
        vkey.valid = true;

        const bool value_changed = !(m_cache.vkey == vkey);
        if (value_changed)
        {
          m_cache.vkey = vkey;

          const auto& rq = qf.getPoint(qp);

          // Cheap scalar P1 element (no allocations)
          const P1Element<ScalarType> fe_scalar(geom);
          const size_t nv = fe_scalar.getCount();

          for (size_t a = 0; a < nv; ++a)
            m_cache.phi_vertex[a] = fe_scalar.getBasis(a)(rq);

          if constexpr (std::is_same_v<RangeType, ScalarType>)
          {
            for (size_t a = 0; a < nv; ++a)
              m_cache.basis[a] = m_cache.phi_vertex[a];
          }
          else
          {
            static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);

            // Update only the one active component; others remain zero.
            for (size_t a = 0; a < nv; ++a)
            {
              const ScalarType val = m_cache.phi_vertex[a];
              for (size_t c = 0; c < vdim; ++c)
              {
                const size_t local = a * vdim + c;
                auto& b = m_cache.basis[local];
                b.coeffRef(c) = val;
              }
            }
          }
        }

        return *this;
      }

      constexpr
      const RangeType& getBasis(size_t local) const
      {
        assert(m_cache.skey);
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

