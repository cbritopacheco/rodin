/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ShapeFunction.h
 * @brief Shape function specializations for the P0 finite element space.
 */
#ifndef RODIN_VARIATIONAL_P0_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_P0_SHAPEFUNCTION_H

#include <utility>
#include <vector>

#include "Rodin/Variational/P0/ForwardDecls.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Math/Traits.h"

namespace Rodin::Variational
{
  /// @brief Shape function expression.
  template <class Derived, class Range, class Mesh, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, P0<Range, Mesh>, Space>
    : public ShapeFunctionBase<ShapeFunction<Derived, P0<Range, Mesh>, Space>, P0<Range, Mesh>, Space>
  {
    public:
      /// @brief Finite element space type.
      using FESType = P0<Range, Mesh>;
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
          /// @brief Vector dimension of the finite element space.
              size_t vdim = 1;
          /// @brief Whether the key holds a cached entry.
              bool valid = false;

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
                return geom == o.geom && vdim == o.vdim;
              }

          /**
           * @brief Resets the key, invalidating the cached entry.
           * @param other Object to copy from.
           */
              void operator=(std::initializer_list<int> other) noexcept
              {
                valid = false;
                geom = Geometry::Polytope::Type::Point;
                vdim = 1;
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
      constexpr
      ShapeFunction(const FESType& fes)
        : Parent(fes),
          m_ip(nullptr)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      constexpr
      ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_ip(nullptr),
          m_cache(other.m_cache)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      constexpr
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_cache(std::move(other.m_cache))
      {}

      /**
       * @brief Returns the number of local basis functions for a polytope.
       * @param polytope Mesh entity used by this operation.
       * @returns Number of local basis functions on the selected entity.
       */
      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        if constexpr (std::is_same_v<RangeType, ScalarType>)
        {
          return P0Element<ScalarType>(polytope.getGeometry()).getCount();
        }
        else
        {
          static_assert(FormLanguage::IsVectorRange<RangeType>::Value);
          return P0Element<RangeType>(polytope.getGeometry(), this->getFiniteElementSpace().getVectorDimension()).getCount();
        }
      }

      /**
       * @brief Gets the integration point the expression is evaluated at.
       * @returns The integration point the expression is evaluated at.
       */
      constexpr
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /**
       * @brief Sets the integration point the expression is evaluated at.
       * @param ip Integration point at which the expression is evaluated.
       * @returns Reference to this object after the operation.
       */
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
            static_assert(FormLanguage::IsVectorRange<RangeType>::Value);

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

      /**
       * @brief Gets the basis function of a local degree of freedom.
       * @param local Index in the local numbering.
       * @returns Value of the selected local basis function at the evaluation point.
       */
      constexpr
      const RangeType& getBasis(size_t local) const
      {
        assert(m_cache.key);
        assert(local < m_cache.basis.size());
        return m_cache.basis[local];
      }

      /**
       * @brief Gets the operand in the shape function expression.
       * @returns The operand in the shape function expression.
       */
      constexpr
      const auto& getLeaf() const
      {
        return static_cast<const Derived&>(*this).getLeaf();
      }

      /**
       * @brief Returns the polynomial order used on a mesh entity.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       * @param polytope Mesh entity; the reported order is independent of this argument.
       */
      constexpr Optional<size_t> getOrder(
        [[maybe_unused]] const Geometry::Polytope& polytope) const noexcept
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

namespace Rodin::Variational
{
  /// @brief Matrix-valued shape-function specialization with cached scalar basis tabulation.
  template <class Derived, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, P0<Math::SpatialMatrix<Scalar>, Mesh>, Space>
    : public ShapeFunctionBase<
        ShapeFunction<Derived, P0<Math::SpatialMatrix<Scalar>, Mesh>, Space>,
        P0<Math::SpatialMatrix<Scalar>, Mesh>, Space>
  {
    public:
      /// @brief CRTP or finite element base class.
      using FES = P0<Math::SpatialMatrix<Scalar>, Mesh>;
      /// @brief Current shape-function specialization.
      using Shape = ShapeFunction;
      /// @brief Existing shape-function interface.
      using Parent = ShapeFunctionBase<
        ShapeFunction<Derived, P0<Math::SpatialMatrix<Scalar>, Mesh>, Space>, FES, Space>;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = typename FormLanguage::Traits<FES>::RangeType;
      /**
       * @brief Constructs matrix basis tabulation on the supplied finite element space.
       * @param fes Finite element space.
       */
      explicit ShapeFunction(const FES& fes)
        : Parent(fes)
      {}
      /**
       * @brief Constructs matrix basis tabulation on the supplied finite element space.
       * @param other Object to copy from.
       */
      ShapeFunction(const ShapeFunction& other)
        : Parent(other)
      {}
      /**
       * @brief Constructs matrix basis tabulation on the supplied finite element space.
       * @param other Object to move from.
       */
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other))
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
       * @brief Binds the integration point and prepares local basis values.
       * @param ip Integration point at which the expression is evaluated.
       * @returns Reference to this object after the operation.
       */
      Shape& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;
        const auto& poly = ip.getPoint().getPolytope();
        const auto& fe = this->getFiniteElementSpace().getFiniteElement(
          poly.getDimension(), poly.getIndex());
        if (!ip.getQuadratureFormula() || m_qf != ip.getQuadratureFormula() ||
          m_qp != ip.getIndex() || m_geometry != poly.getGeometry())
        {
          m_basis.resize(fe.getCount());
          const auto& scalar = fe.getScalarElement();
          const size_t components = fe.getRows() * fe.getColumns();
          const auto& point = ip.getPoint().getReferenceCoordinates();
          for (size_t a = 0; a < scalar.getCount(); ++a)
          {
            const auto value = scalar.getBasis(a)(point);
            for (size_t c = 0; c < components; ++c)
            {
              auto& basis = m_basis[a * components + c];
              basis.resize(fe.getRows(), fe.getColumns());
              basis.setZero();
              basis(c / fe.getColumns(), c % fe.getColumns()) = value;
            }
          }
          m_qf = ip.getQuadratureFormula();
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
        return this->getFiniteElementSpace()
          .getFiniteElement(poly.getDimension(), poly.getIndex())
          .getOrder();
      }

      ShapeFunction* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      const IntegrationPoint* m_ip = nullptr;
      const QF::QuadratureFormulaBase* m_qf = nullptr;
      size_t m_qp = 0;
      Geometry::Polytope::Type m_geometry = Geometry::Polytope::Type::Point;
      std::vector<RangeType> m_basis;
  };
}

#endif
