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
    : public ShapeFunctionBase<
        ShapeFunction<Derived, H1<K, Scalar, Mesh>, Space>,
        H1<K, Scalar, Mesh>,
        Space>
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

      struct Cache
      {
        struct Key
        {
          Geometry::Polytope::Type geom = Geometry::Polytope::Type::Point;
          const QF::QuadratureFormulaBase* qf = nullptr;
          size_t qp = 0;
          bool valid = false;

          explicit operator bool() const noexcept { return valid; }

          bool operator==(const Key& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return geom == o.geom && qf == o.qf && qp == o.qp;
          }

          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            geom = Geometry::Polytope::Type::Point;
            qf = nullptr;
            qp = 0;
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
          m_ip(nullptr)
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
        // Fast, geometry-only.
        return H1Element<K, ScalarType>(polytope.getGeometry()).getCount();
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

        const auto& p   = ip.getPoint();
        const auto& poly = p.getPolytope();
        const auto  geom = poly.getGeometry();

        const auto& qf  = ip.getQuadratureFormula();
        const size_t qp = ip.getIndex();

        typename Cache::Key key;
        key.geom  = geom;
        key.qf    = &qf;
        key.qp    = qp;
        key.valid = true;

        if (!(m_cache.key == key))
        {
          m_cache.key = key;

          const H1Element<K, ScalarType> fe(geom);
          const size_t ndof = fe.getCount();

          m_cache.basis.resize(ndof);

          // Fast path: use element tabulation (basis + gradients in ref coords)
          const auto& tab = fe.getTabulation(qf);

          for (size_t a = 0; a < ndof; ++a)
            m_cache.basis[a] = tab.getBasis(qp, a);
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
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        // Uses your H1Element::getOrder() (K / 2K / 3K depending on geometry)
        return H1Element<K, ScalarType>(geom.getGeometry()).getOrder();
      }

      ShapeFunction* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      const IntegrationPoint* m_ip;
      Cache m_cache;
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

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType; // == Scalar
      using RangeType  = typename FormLanguage::Traits<FESType>::RangeType;  // == Math::Vector<Scalar>

      using Parent =
        ShapeFunctionBase<
          ShapeFunction<Derived, FESType, SpaceType>,
          FESType,
          SpaceType>;

      static_assert(std::is_same_v<RangeType, Math::Vector<ScalarType>>);

      struct Cache
      {
        struct Key
        {
          Geometry::Polytope::Type geom = Geometry::Polytope::Type::Point;
          const QF::QuadratureFormulaBase* qf = nullptr;
          size_t qp = 0;
          size_t vdim = 0;
          bool valid = false;

          explicit operator bool() const noexcept { return valid; }

          bool operator==(const Key& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return geom == o.geom && qf == o.qf && qp == o.qp && vdim == o.vdim;
          }

          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            geom = Geometry::Polytope::Type::Point;
            qf = nullptr;
            qp = 0;
            vdim = 0;
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
        const size_t vdim = this->getFiniteElementSpace().getVectorDimension();
        const size_t ndof_scalar = H1Element<K, ScalarType>(polytope.getGeometry()).getCount();
        return ndof_scalar * vdim;
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

        const auto& p    = ip.getPoint();
        const auto& poly = p.getPolytope();
        const auto  geom = poly.getGeometry();

        const auto& qf  = ip.getQuadratureFormula();
        const size_t qp = ip.getIndex();

        const size_t vdim = this->getFiniteElementSpace().getVectorDimension();

        typename Cache::Key key;
        key.geom  = geom;
        key.qf    = &qf;
        key.qp    = qp;
        key.vdim  = vdim;
        key.valid = true;

        if (!(m_cache.key == key))
        {
          m_cache.key = key;

          // Scalar tabulation
          const H1Element<K, ScalarType> fe_scalar(geom);
          const size_t ndof_scalar = fe_scalar.getCount();
          const size_t ndof = ndof_scalar * vdim;

          m_cache.basis.resize(ndof);

          const auto& tab = fe_scalar.getTabulation(qf);

          // φ_{a,c} = φ_a e_c
          for (size_t a = 0; a < ndof_scalar; ++a)
          {
            const ScalarType val = tab.getBasis(qp, a);
            for (size_t c = 0; c < vdim; ++c)
            {
              RangeType v;
              v.resize(vdim);
              v.setZero();
              v.coeffRef(c) = val;
              m_cache.basis[a * vdim + c] = std::move(v);
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
      Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        // Order of the underlying scalar polynomial space (your total-degree convention)
        return H1Element<K, ScalarType>(geom.getGeometry()).getOrder();
      }

      ShapeFunction* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      const IntegrationPoint* m_ip;
      Cache m_cache;
  };
}

#endif
