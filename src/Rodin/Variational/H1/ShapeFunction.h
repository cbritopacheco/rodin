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
  template <class Derived, size_t K, class Scalar, class Mesh,
    ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, H1<K, Scalar, Mesh>, Space>
    : public ShapeFunctionBase<ShapeFunction<Derived, H1<K, Scalar, Mesh>, Space>,
        H1<K, Scalar, Mesh>, Space>
  {
    public:
      /// @brief Finite element space type.
      using FESType = H1<K, Scalar, Mesh>;

      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Range (evaluation value) type.
      using RangeType = typename FormLanguage::Traits<FESType>::RangeType;

      /// @brief Parent class type.
      using Parent =
        ShapeFunctionBase<ShapeFunction<Derived, FESType, SpaceType>, FESType, SpaceType>;

      /// @brief Per-cell tabulation cache.
      struct Cache
      {
        /// @brief Key identifying a cached tabulation.
        struct Key
        {
          /// @brief Geometry of the cached polytope.
          Geometry::Polytope::Type geom;

          /// @brief Quadrature formula the cached tabulation belongs to.
          const QF::QuadratureFormula* qf;

          /// @brief Index of the quadrature point.
          size_t qp;

          /// @brief Whether the key holds a cached entry.
          bool valid;

          /// @brief Constructs an invalid quadrature cache key.
          Key()
            : geom(Geometry::Polytope::Type::Point),
              qf(nullptr),
              qp(0),
              valid(false)
          {}

          /**
           * @brief Tests whether the key holds a cached entry.
           * @returns True if the key identifies a cached entry; false otherwise.
           */
          explicit operator bool() const noexcept
          {
            return valid;
          }

          /**
           * @brief Equality comparison.
           * @returns Whether the operands compare equal.
           * @param o Key to compare with this key.
           */
          bool operator==(const Key& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return geom == o.geom && qf == o.qf && qp == o.qp;
          }

          /// @brief Resets the key, invalidating the cached entry.
          void reset() noexcept
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

      /**
       * @brief Constructs the shape function over a finite element space.
       * @param fes Finite element space.
       */
      constexpr ShapeFunction(const FESType& fes)
        : Parent(fes),
          m_ip(nullptr)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      constexpr ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_ip(nullptr)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      constexpr ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_cache(std::move(other.m_cache))
      {}

      /**
       * @brief Returns the number of local basis functions for a polytope.
       * @param polytope Mesh entity used by this operation.
       * @returns Number of local basis functions on the selected entity.
       */
      constexpr size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        // Fast, geometry-only.
        return H1Element<K, ScalarType>(polytope.getGeometry()).getCount();
      }

      /**
       * @brief Gets the integration point the expression is evaluated at.
       * @returns The integration point the expression is evaluated at.
       */
      constexpr const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /**
       * @brief Starts a fresh cell binding and releases the previous point binding.
       * @param polytope Entity whose basis will be evaluated.
       * @returns This expression after rebinding its evaluation state.
       */
      ShapeFunction& setPolytope([[maybe_unused]] const Geometry::Polytope& polytope)
      {
        m_ip = nullptr;
        return *this;
      }

      /**
       * @brief Sets the integration point the expression is evaluated at.
       *
       * Pointwise interpolation-node bindings use exact cardinal values.
       * This preserves structural zeros in identification constraints; nearby
       * points and quadrature samples retain ordinary basis evaluation.
       * @param ip Integration point at which the expression is evaluated.
       * @returns Reference to this object after the operation.
       */
      ShapeFunction& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& p = ip.getPoint();
        const auto& poly = p.getPolytope();
        const auto geom = poly.getGeometry();

        const auto* qf = ip.getReferenceQuadratureFormula();
        const size_t qp = qf ? ip.getIndex() : 0;

        typename Cache::Key key;
        key.geom = geom;
        key.qf = qf;
        key.qp = qp;
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
            const auto node = fe.getNodeIndex(rc);
            for (size_t a = 0; a < ndof; ++a)
            {
              m_cache.basis[a] = node ? ScalarType(a == *node) : fe.getBasis(a)(rc);
            }
          }
        }

        return *this;
      }

      /**
       * @brief Gets the basis function of a local degree of freedom.
       * @param local Index in the local numbering.
       * @returns Value of the selected local basis function at the evaluation point.
       */
      constexpr const RangeType& getBasis(size_t local) const
      {
        assert(m_cache.key);
        assert(local < m_cache.basis.size());
        return m_cache.basis[local];
      }

      /**
       * @brief Gets the operand in the shape function expression.
       * @returns The operand in the shape function expression.
       */
      constexpr const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      /**
       * @brief Returns the polynomial order used on a mesh entity.
       * @param geom Reference geometry.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      constexpr Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        // Retain the total-degree estimate; rational pyramids use a family accuracy policy.
        return H1Element<K, ScalarType>(geom.getGeometry()).getDegree().value_or(2 * K);
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
  template <class Derived, size_t K, class Scalar, class Mesh,
    ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, H1<K, Math::SpatialVector<Scalar>, Mesh>, Space>
    : public ShapeFunctionBase<
        ShapeFunction<Derived, H1<K, Math::SpatialVector<Scalar>, Mesh>, Space>,
        H1<K, Math::SpatialVector<Scalar>, Mesh>, Space>
  {
    public:
      /// @brief Finite element space type.
      using FESType = H1<K, Math::SpatialVector<Scalar>, Mesh>;

      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType; // == Scalar

      /// @brief Range (evaluation value) type.
      using RangeType =
        typename FormLanguage::Traits<FESType>::RangeType; // == Math::Vector<Scalar>

      /// @brief Parent class type.
      using Parent =
        ShapeFunctionBase<ShapeFunction<Derived, FESType, SpaceType>, FESType, SpaceType>;

      static_assert(FormLanguage::IsVectorRange<RangeType>::Value);

      /// @brief Per-cell tabulation cache.
      struct Cache
      {
        /// @brief Key identifying a cached tabulation.
        struct Key
        {
          /// @brief Geometry of the cached polytope.
          Geometry::Polytope::Type geom;

          /// @brief Quadrature formula the cached tabulation belongs to.
          const QF::QuadratureFormula* qf;

          /// @brief Index of the quadrature point.
          size_t qp;

          /// @brief Vector dimension of the finite element space.
          size_t vdim;

          /// @brief Whether the key holds a cached entry.
          bool valid;

          /// @brief Constructs an invalid quadrature cache key.
          Key()
            : geom(Geometry::Polytope::Type::Point),
              qf(nullptr),
              qp(0),
              vdim(0),
              valid(false)
          {}

          /**
           * @brief Tests whether the key holds a cached entry.
           * @returns True if the key identifies a cached entry; false otherwise.
           */
          explicit operator bool() const noexcept
          {
            return valid;
          }

          /**
           * @brief Equality comparison.
           * @returns Whether the operands compare equal.
           * @param o Key to compare with this key.
           */
          bool operator==(const Key& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return geom == o.geom && qf == o.qf && qp == o.qp && vdim == o.vdim;
          }

          /// @brief Resets the key, invalidating the cached entry.
          void reset() noexcept
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

      /**
       * @brief Constructs the shape function over a finite element space.
       * @param fes Finite element space.
       */
      constexpr ShapeFunction(const FESType& fes)
        : Parent(fes),
          m_ip(nullptr)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      constexpr ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_ip(nullptr),
          m_cache(other.m_cache)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      constexpr ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_cache(std::move(other.m_cache))
      {}

      /**
       * @brief Returns the number of local basis functions for a polytope.
       * @param polytope Mesh entity used by this operation.
       * @returns Number of local basis functions on the selected entity.
       */
      constexpr size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        const size_t vdim = this->getFiniteElementSpace().getVectorDimension();
        const size_t ndofScalar =
          H1Element<K, ScalarType>(polytope.getGeometry()).getCount();
        return ndofScalar * vdim;
      }

      /**
       * @brief Gets the integration point the expression is evaluated at.
       * @returns The integration point the expression is evaluated at.
       */
      constexpr const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /**
       * @brief Starts a fresh cell binding and releases the previous point binding.
       * @param polytope Entity whose basis will be evaluated.
       * @returns This expression after rebinding its evaluation state.
       */
      ShapeFunction& setPolytope([[maybe_unused]] const Geometry::Polytope& polytope)
      {
        m_ip = nullptr;
        return *this;
      }

      /**
       * @brief Sets the integration point the expression is evaluated at.
       * @param ip Integration point at which the expression is evaluated.
       * @returns Reference to this object after the operation.
       */
      ShapeFunction& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& p = ip.getPoint();
        const auto& poly = p.getPolytope();
        const auto geom = poly.getGeometry();

        const auto* qf = ip.getReferenceQuadratureFormula();
        const size_t qp = qf ? ip.getIndex() : 0;

        const size_t vdim = this->getFiniteElementSpace().getVectorDimension();

        typename Cache::Key key;
        key.geom = geom;
        key.qf = qf;
        key.qp = qp;
        key.vdim = vdim;
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
          const auto node = qf ? Optional<size_t>{} : feScalar.getNodeIndex(rc);

          // φ_{a,c} = φ_a e_c
          for (size_t a = 0; a < ndofScalar; ++a)
          {
            const ScalarType val = qf ? tab->getBasis(qp, a)
              : node                  ? ScalarType(a == *node)
                                      : feScalar.getBasis(a)(rc);
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

      /**
       * @brief Gets the basis function of a local degree of freedom.
       * @param local Index in the local numbering.
       * @returns Value of the selected local basis function at the evaluation point.
       */
      constexpr const RangeType& getBasis(size_t local) const
      {
        assert(m_cache.key);
        assert(local < m_cache.basis.size());
        return m_cache.basis[local];
      }

      /**
       * @brief Gets the operand in the shape function expression.
       * @returns The operand in the shape function expression.
       */
      constexpr const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      /**
       * @brief Returns the polynomial order used on a mesh entity.
       * @param geom Reference geometry.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      constexpr Optional<size_t> getOrder(const Geometry::Polytope& geom) const noexcept
      {
        // Order of the underlying scalar polynomial space (your total-degree convention)
        return H1Element<K, ScalarType>(geom.getGeometry()).getDegree().value_or(2 * K);
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

namespace Rodin::Variational
{
  /// @brief Matrix-valued shape-function specialization with cached scalar basis tabulation.
  template <class Derived, size_t K, class Scalar, class Mesh,
    ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, H1<K, Math::SpatialMatrix<Scalar>, Mesh>, Space>
    : public ShapeFunctionBase<
        ShapeFunction<Derived, H1<K, Math::SpatialMatrix<Scalar>, Mesh>, Space>,
        H1<K, Math::SpatialMatrix<Scalar>, Mesh>, Space>
  {
    public:
      /// @brief CRTP or finite element base class.
      using FES = H1<K, Math::SpatialMatrix<Scalar>, Mesh>;

      /// @brief Current shape-function specialization.
      using Shape = ShapeFunction;

      /// @brief Existing shape-function interface.
      using Parent = ShapeFunctionBase<
        ShapeFunction<Derived, H1<K, Math::SpatialMatrix<Scalar>, Mesh>, Space>, FES,
        Space>;

      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = typename FormLanguage::Traits<FES>::RangeType;

      /**
       * @brief Constructs matrix basis tabulation on the supplied finite element space.
       * @param fes Finite element space.
       */
      explicit ShapeFunction(const FES& fes)
        : Parent(fes),
          m_ip(nullptr),
          m_qf(nullptr),
          m_qp(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Constructs matrix basis tabulation on the supplied finite element space.
       * @param other Object to copy from.
       */
      ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_ip(nullptr),
          m_qf(nullptr),
          m_qp(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Constructs matrix basis tabulation on the supplied finite element space.
       * @param other Object to move from.
       */
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_ip(nullptr),
          m_qf(nullptr),
          m_qp(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      /**
       * @brief Returns the local basis count for the selected polytope.
       * @param poly Mesh entity used by this operation.
       * @returns Number of local basis functions on the selected entity.
       */
      size_t getDOFs(const Geometry::Polytope& poly) const
      {
        return this->getFiniteElementSpace()
          .getFiniteElement(poly.getDimension(), poly.getIndex())
          .getCount();
      }

      /**
       * @brief Returns the currently bound integration point.
       * @returns The currently bound integration point.
       */
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /**
       * @brief Starts a fresh cell binding and releases the previous point binding.
       * @param polytope Entity whose basis will be evaluated.
       * @returns This expression after rebinding its evaluation state.
       */
      Shape& setPolytope([[maybe_unused]] const Geometry::Polytope& polytope)
      {
        m_ip = nullptr;
        return *this;
      }

      /**
       * @brief Binds the integration point and prepares local basis values.
       *
       * Cached component values depend on both matrix dimensions, not only
       * their product. Reassigning a space can change @f$2\times3@f$ to
       * @f$3\times2@f$ without changing the local basis count or sample identity.
       * The owned basis layout therefore participates in cache invalidation.
       *
       * A certified reference sample may borrow an already prepared scalar
       * table during this call. A missing or evicted table retains direct
       * pointwise evaluation: requesting one matrix value must not construct
       * values and derivatives for an entire formula. No table borrow is
       * retained after the owned component basis has been filled.
       *
       * @param ip Integration point at which the expression is evaluated.
       * @returns Reference to this object after the operation.
       */
      Shape& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;
        const auto& poly = ip.getPoint().getPolytope();
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(
          poly.getDimension(), poly.getIndex());
        if (!ip.getReferenceQuadratureFormula() ||
          m_qf != ip.getReferenceQuadratureFormula() || m_qp != ip.getIndex() ||
          m_geometry != poly.getGeometry() || m_basis.size() != fe.getCount() ||
          m_basis.empty() || m_basis.front().rows() != fe.getRows() ||
          m_basis.front().cols() != fe.getColumns())
        {
          m_basis.resize(fe.getCount());
          const auto& scalar = fe.getScalarElement();
          const size_t components = fe.getRows() * fe.getColumns();
          const auto& point = ip.getPoint().getReferenceCoordinates();
          Optional<
            std::reference_wrapper<const typename H1Element<K, Scalar>::Tabulation>>
            tab;
          if (const auto* formula = ip.getReferenceQuadratureFormula())
            tab = scalar.findTabulation(*formula);
          const auto node = ip.getReferenceQuadratureFormula()
            ? Optional<size_t>{}
            : scalar.getNodeIndex(point);
          for (size_t a = 0; a < scalar.getCount(); ++a)
          {
            const auto value = tab ? tab->get().getBasis(ip.getIndex(), a)
              : node               ? Scalar(a == *node)
                                   : scalar.getBasis(a)(point);
            for (size_t c = 0; c < components; ++c)
            {
              auto& basis = m_basis[a * components + c];
              basis.resize(fe.getRows(), fe.getColumns());
              basis.setZero();
              basis(c / fe.getColumns(), c % fe.getColumns()) = value;
            }
          }
          m_qf = ip.getReferenceQuadratureFormula();
          m_qp = ip.getIndex();
          m_geometry = poly.getGeometry();
        }
        return static_cast<Shape&>(*this);
      }

      /**
       * @brief Returns a basis value at the bound integration point.
       * @param local Index in the local numbering.
       * @returns Value of the selected local basis function at the evaluation point.
       */
      const RangeType& getBasis(size_t local) const
      {
        return m_basis.at(local);
      }

      /**
       * @brief Returns the leaf shape function used for assembly.
       * @returns The leaf shape function used for assembly.
       */
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      /**
       * @brief Returns the polynomial order when it is known.
       * @param poly Mesh entity used by this operation.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        const auto& element = this->getFiniteElementSpace().getFiniteElement(
          poly.getDimension(), poly.getIndex());
        return element.getDegree().value_or(2 * K);
      }

      ShapeFunction* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      const IntegrationPoint* m_ip;
      const QF::QuadratureFormula* m_qf;
      size_t m_qp;
      Geometry::Polytope::Type m_geometry;
      std::vector<RangeType> m_basis;
  };
}

#endif
