/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ShapeFunction.h
 * @brief Shape function specializations for the H1 finite element space.
 */
#ifndef RODIN_VARIATIONAL_H1_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_H1_SHAPEFUNCTION_H

#include <utility>
#include <vector>

#include "Rodin/Variational/H1/ForwardDecls.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Math/Traits.h"

namespace Rodin::Variational
{
  /// @brief Shape function expression.
  template <class Derived, size_t K, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, H1<K, Scalar, Mesh>, Space>
    : public ShapeFunctionBase<
        ShapeFunction<Derived, H1<K, Scalar, Mesh>, Space>,
        H1<K, Scalar, Mesh>,
        Space>
  {
    public:
      /// @brief Finite element space type.
      using FESType = H1<K, Scalar, Mesh>;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      /// @brief Range (evaluation value) type.
      using RangeType  = typename FormLanguage::Traits<FESType>::RangeType;

      /// @brief Parent class type.
      using Parent =
        ShapeFunctionBase<
          ShapeFunction<Derived, FESType, SpaceType>,
          FESType,
          SpaceType>;

      /// @brief Per-cell tabulation cache.
      struct Cache
      {
        /// @brief Key identifying a cached tabulation.
        struct Key
        {
          /// @brief Geometry of the cached polytope.
          Geometry::Polytope::Type geom = Geometry::Polytope::Type::Point;
          /// @brief Quadrature formula the cached tabulation belongs to.
          const QF::QuadratureFormulaBase* qf = nullptr;
          /// @brief Index of the quadrature point.
          size_t qp = 0;
          /// @brief Whether the key holds a cached entry.
          bool valid = false;

          /// @brief Tests whether the key holds a cached entry.
          explicit operator bool() const noexcept { return valid; }

          /// @brief Equality comparison.
          bool operator==(const Key& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return geom == o.geom && qf == o.qf && qp == o.qp;
          }

          /// @brief Resets the key, invalidating the cached entry.
          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            geom = Geometry::Polytope::Type::Point;
            qf = nullptr;
            qp = 0;
          }
        };

        /// @brief Cached reference basis tabulation.
        std::vector<RangeType> basis;
        /// @brief Key identifying the cached entry.
        Key key;
      };

      ShapeFunction() = delete;

      constexpr
      /// @brief Constructs the shape function over a finite element space.
      ShapeFunction(const FESType& fes)
        : Parent(fes),
          m_ip(nullptr)
      {}

      constexpr
      /// @brief Copy constructor.
      ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_ip(nullptr)
      {}

      constexpr
      /// @brief Move constructor.
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_cache(std::move(other.m_cache))
      {}

      constexpr
      /// @brief Gets the global DOF indices for a polytope.
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        // Fast, geometry-only.
        return H1Element<K, ScalarType>(polytope.getGeometry()).getCount();
      }

      constexpr
      /// @brief Gets the integration point the expression is evaluated at.
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /// @brief Sets the integration point the expression is evaluated at.
      ShapeFunction& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& p   = ip.getPoint();
        const auto& poly = p.getPolytope();
        const auto  geom = poly.getGeometry();

        const auto* qf = ip.getQuadratureFormula();
        const size_t qp = qf ? ip.getIndex() : 0;

        typename Cache::Key key;
        key.geom  = geom;
        key.qf = qf;
        key.qp    = qp;
        key.valid = true;

        if (!qf || !(m_cache.key == key))
        {
          m_cache.key = key;

          const H1Element<K, ScalarType> fe(geom);
          const size_t ndof = fe.getCount();

          m_cache.basis.resize(ndof);

          if (qf)
          {
            // Fast path: use element tabulation inside quadrature loops.
            const auto& tab = fe.getTabulation(*qf);
            for (size_t a = 0; a < ndof; ++a)
              m_cache.basis[a] = tab.getBasis(qp, a);
          }
          else
          {
            const auto& rc = p.getReferenceCoordinates();
            for (size_t a = 0; a < ndof; ++a)
              m_cache.basis[a] = fe.getBasis(a)(rc);
          }
        }

        return *this;
      }

      constexpr
      /// @brief Gets the basis function of a local degree of freedom.
      const RangeType& getBasis(size_t local) const
      {
        assert(m_cache.key);
        assert(local < m_cache.basis.size());
        return m_cache.basis[local];
      }

      constexpr
      /// @brief Gets the operand in the shape function expression.
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      constexpr
      /// @brief Returns the polynomial order used on a mesh entity.
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

  /// @brief Shape function of a vector-valued H1 finite element space.
  template <class Derived, size_t K, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, H1<K, Math::SpatialVector<Scalar>, Mesh>, Space>
    : public ShapeFunctionBase<
        ShapeFunction<Derived, H1<K, Math::SpatialVector<Scalar>, Mesh>, Space>,
        H1<K, Math::SpatialVector<Scalar>, Mesh>,
        Space>
  {
    public:
      /// @brief Finite element space type.
      using FESType = H1<K, Math::SpatialVector<Scalar>, Mesh>;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType; // == Scalar
      /// @brief Range (evaluation value) type.
      using RangeType  = typename FormLanguage::Traits<FESType>::RangeType;  // == Math::Vector<Scalar>

      /// @brief Parent class type.
      using Parent =
        ShapeFunctionBase<
          ShapeFunction<Derived, FESType, SpaceType>,
          FESType,
          SpaceType>;

      static_assert(FormLanguage::IsVectorRange<RangeType>::Value);

      /// @brief Per-cell tabulation cache.
      struct Cache
      {
        /// @brief Key identifying a cached tabulation.
        struct Key
        {
          /// @brief Geometry of the cached polytope.
          Geometry::Polytope::Type geom = Geometry::Polytope::Type::Point;
          /// @brief Quadrature formula the cached tabulation belongs to.
          const QF::QuadratureFormulaBase* qf = nullptr;
          /// @brief Index of the quadrature point.
          size_t qp = 0;
          /// @brief Vector dimension of the finite element space.
          size_t vdim = 0;
          /// @brief Whether the key holds a cached entry.
          bool valid = false;

          /// @brief Tests whether the key holds a cached entry.
          explicit operator bool() const noexcept { return valid; }

          /// @brief Equality comparison.
          bool operator==(const Key& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return geom == o.geom && qf == o.qf && qp == o.qp && vdim == o.vdim;
          }

          /// @brief Resets the key, invalidating the cached entry.
          void operator=(std::initializer_list<int>) noexcept
          {
            valid = false;
            geom = Geometry::Polytope::Type::Point;
            qf = nullptr;
            qp = 0;
            vdim = 0;
          }
        };

        /// @brief Cached reference basis tabulation.
        std::vector<RangeType> basis;
        /// @brief Key identifying the cached entry.
        Key key;
      };

      ShapeFunction() = delete;

      constexpr
      /// @brief Constructs the shape function over a finite element space.
      ShapeFunction(const FESType& fes)
        : Parent(fes),
          m_ip(nullptr)
      {}

      constexpr
      /// @brief Copy constructor.
      ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_ip(nullptr),
          m_cache(other.m_cache)
      {}

      constexpr
      /// @brief Move constructor.
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_cache(std::move(other.m_cache))
      {}

      constexpr
      /// @brief Gets the global DOF indices for a polytope.
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        const size_t vdim = this->getFiniteElementSpace().getVectorDimension();
        const size_t ndofScalar =
          H1Element<K, ScalarType>(polytope.getGeometry()).getCount();
        return ndofScalar * vdim;
      }

      constexpr
      /// @brief Gets the integration point the expression is evaluated at.
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /// @brief Sets the integration point the expression is evaluated at.
      ShapeFunction& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& p    = ip.getPoint();
        const auto& poly = p.getPolytope();
        const auto  geom = poly.getGeometry();

        const auto* qf = ip.getQuadratureFormula();
        const size_t qp = qf ? ip.getIndex() : 0;

        const size_t vdim = this->getFiniteElementSpace().getVectorDimension();

        typename Cache::Key key;
        key.geom  = geom;
        key.qf = qf;
        key.qp    = qp;
        key.vdim  = vdim;
        key.valid = true;

        if (!qf || !(m_cache.key == key))
        {
          m_cache.key = key;

          // Scalar tabulation
          const H1Element<K, ScalarType> feScalar(geom);
          const size_t ndofScalar = feScalar.getCount();
          const size_t ndof = ndofScalar * vdim;

          m_cache.basis.resize(ndof);

          const auto* tab = qf ? &feScalar.getTabulation(*qf) : nullptr;
          const auto& rc = p.getReferenceCoordinates();

          // φ_{a,c} = φ_a e_c
          for (size_t a = 0; a < ndofScalar; ++a)
          {
            const ScalarType val = qf ? tab->getBasis(qp, a) : feScalar.getBasis(a)(rc);
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
      /// @brief Gets the basis function of a local degree of freedom.
      const RangeType& getBasis(size_t local) const
      {
        assert(m_cache.key);
        assert(local < m_cache.basis.size());
        return m_cache.basis[local];
      }

      constexpr
      /// @brief Gets the operand in the shape function expression.
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      constexpr
      /// @brief Returns the polynomial order used on a mesh entity.
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
