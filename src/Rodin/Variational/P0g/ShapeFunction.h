#ifndef RODIN_VARIATIONAL_P0G_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_P0G_SHAPEFUNCTION_H

#include <type_traits>
#include <utility>
#include <vector>

#include "Rodin/Variational/P0g/ForwardDecls.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/IntegrationPoint.h"

namespace Rodin::Variational
{
  template <class Derived, class Range, class Mesh, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, P0g<Range, Mesh>, Space>
    : public ShapeFunctionBase<
        ShapeFunction<Derived, P0g<Range, Mesh>, Space>,
        P0g<Range, Mesh>,
        Space>
  {
    public:
      using FESType = P0g<Range, Mesh>;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      using RangeType  = typename FormLanguage::Traits<FESType>::RangeType;

      using Parent =
        ShapeFunctionBase<
          ShapeFunction<Derived, FESType, SpaceType>,
          FESType,
          SpaceType>;

      static constexpr bool IsScalarRange = std::is_same_v<RangeType, ScalarType>;
      static constexpr bool IsVectorRange = std::is_same_v<RangeType, Math::Vector<ScalarType>>;

      static_assert(IsScalarRange || IsVectorRange);

      // Cache is only instantiated/used for the vector case.
      struct VectorCache
      {
        size_t vdim = 0;
        std::vector<RangeType> basis; // basis[c] = e_c (size vdim)
      };

      ShapeFunction() = delete;

      constexpr
      explicit ShapeFunction(const FESType& fes)
        : Parent(fes),
          m_ip(nullptr)
      {}

      constexpr
      ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_ip(nullptr),
          m_vcache(other.m_vcache)
      {}

      constexpr
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_vcache(std::move(other.m_vcache))
      {}

      constexpr
      size_t getDOFs(const Geometry::Polytope&) const
      {
        if constexpr (IsScalarRange)
        {
          return 1;
        }
        else
        {
          static_assert(IsVectorRange);
          return this->getFiniteElementSpace().getVectorDimension();
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
        // P0g basis does not depend on the integration point.
        m_ip = &ip;

        if constexpr (IsVectorRange)
        {
          const size_t vdim = this->getFiniteElementSpace().getVectorDimension();
          if (m_vcache.vdim != vdim)
            rebuildVectorCache(vdim);
        }

        return *this;
      }

      constexpr
      decltype(auto) getBasis(size_t local) const
      {
        if constexpr (IsScalarRange)
        {
          assert(local == 0);
          (void) local;
          return ScalarType(1);
        }
        else
        {
          static_assert(IsVectorRange);
          assert(m_vcache.vdim > 0);
          assert(local < m_vcache.basis.size());
          return m_vcache.basis[local];
        }
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      constexpr
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      virtual ShapeFunction* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      void rebuildVectorCache(size_t vdim)
      {
        m_vcache.vdim = vdim;
        m_vcache.basis.resize(vdim);

        for (size_t c = 0; c < vdim; ++c)
        {
          RangeType e;
          e.resize(vdim);
          e.setZero();
          e.coeffRef(c) = ScalarType(1);
          m_vcache.basis[c] = std::move(e);
        }
      }

      const IntegrationPoint* m_ip;

      // Exists for both cases, but only used/rebuilt in the vector case.
      VectorCache m_vcache;
  };
}

#endif
