#ifndef RODIN_VARIATIONAL_P0_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_P0_SHAPEFUNCTION_H

#include <utility>
#include <vector>

#include "Rodin/Variational/P0/ForwardDecls.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/IntegrationPoint.h"

namespace Rodin::Variational
{
  template <class Derived, class Range, class Mesh, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, P0<Range, Mesh>, Space>
    : public ShapeFunctionBase<ShapeFunction<Derived, P0<Range, Mesh>, Space>, P0<Range, Mesh>, Space>
  {
    public:
      using FESType = P0<Range, Mesh>;
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
          size_t vdim = 1;
          bool valid = false;

          explicit operator bool() const noexcept { return valid; }

          bool operator==(const Key& o) const noexcept
          {
            if (!valid || !o.valid) return false;
            return geom == o.geom && vdim == o.vdim;
          }

          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            geom = Geometry::Polytope::Type::Point;
            vdim = 1;
          }
        };

        std::vector<RangeType> basis;
        Key key;
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
          return P0Element<ScalarType>(polytope.getGeometry()).getCount();
        }
        else
        {
          static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);
          return P0Element<RangeType>(polytope.getGeometry(), this->getFiniteElementSpace().getVectorDimension()).getCount();
        }
      }

      constexpr
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      ShapeFunction& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& poly = ip.getPoint().getPolytope();
        const auto geom  = poly.getGeometry();

        const size_t vdim =
          [] (const auto& self) -> size_t
          {
            if constexpr (std::is_same_v<RangeType, ScalarType>)
              return 1;
            else
              return self.getFiniteElementSpace().getVectorDimension();
          }(*this);

        typename Cache::Key key;
        key.geom  = geom;
        key.vdim  = vdim;
        key.valid = true;

        if (!(m_cache.key == key))
        {
          m_cache.key = key;

          if constexpr (std::is_same_v<RangeType, ScalarType>)
          {
            m_cache.basis.resize(1);
            m_cache.basis[0] = ScalarType(1);
          }
          else
          {
            static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);

            m_cache.basis.resize(vdim);
            for (size_t c = 0; c < vdim; ++c)
            {
              RangeType e;
              e.resize(vdim);
              e.setZero();
              e.coeffRef(c) = ScalarType(1);
              m_cache.basis[c] = std::move(e);
            }
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
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        // P0 basis functions are constant on the reference element for any geometry.
        return 0;
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


