/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ShapeFunction.h
 * @brief Shape function specializations for the P1 finite element space.
 */
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
#include "Rodin/Math/Traits.h"

namespace Rodin::Variational
{
  /// @brief Shape function expression.
  template <class Derived, class Range, class Mesh, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, P1<Range, Mesh>, Space>
    : public ShapeFunctionBase<ShapeFunction<Derived, P1<Range, Mesh>, Space>, P1<Range, Mesh>, Space>
  {
    public:
      /// @brief Finite element space type.
      using FESType = P1<Range, Mesh>;
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
        /// @brief Key identifying the cached element structure.
          struct StructureKey
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
              bool operator==(const StructureKey& o) const noexcept
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

        /// @brief Key identifying the cached shape function values.
        struct ValueKey
        {
          /// @brief Quadrature formula the cached tabulation belongs to.
            const QF::QuadratureFormulaBase* qf = nullptr;
          /// @brief Index of the quadrature point.
            size_t qp = 0;
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
            bool operator==(const ValueKey& o) const noexcept
            {
              if (!valid || !o.valid)
                return false;
              return qf == o.qf && qp == o.qp;
            }

          /**
           * @brief Resets the key, invalidating the cached entry.
           * @param other Object to copy from.
           */
            void operator=(std::initializer_list<int> other) noexcept
            {
              valid = false;
              qf = nullptr;
              qp = 0;
            }
        };

        // For scalar: size = nv
        // For vector: size = nv * vdim, where local = a*vdim + c
        /// @brief Cached reference basis tabulation.
        std::vector<RangeType>  basis;

        // Scalar vertex basis values: size = nv
        /// @brief Cached vertex basis values.
        std::vector<ScalarType> phi_vertex;

        /// @brief Key of the cached element structure.
        StructureKey skey;
        /// @brief Key of the cached shape function values.
        ValueKey vkey;
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
          return P1Element<ScalarType>(polytope.getGeometry()).getCount();
        }
        else
        {
          static_assert(FormLanguage::IsVectorRange<RangeType>::Value);
          return P1Element<RangeType>(polytope.getGeometry(), this->getFiniteElementSpace().getVectorDimension()).getCount();
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
       * Fallback path (non-quadrature evaluations).
       * Keeps your previous behavior (pushforward of each basis at p).
       * @param ip Integration point at which the expression is evaluated.
       * @returns Reference to this object after the operation.
       */
      ShapeFunction& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_ip = &ip;

        const auto& p  = ip.getPoint();
        const auto* qf = ip.getQuadratureFormula();
        const size_t qp = qf ? ip.getIndex() : 0;

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

        const bool structureChanged = !(m_cache.skey == skey);
        if (structureChanged)
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
            static_assert(FormLanguage::IsVectorRange<RangeType>::Value);
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
        vkey.qf = qf;
        vkey.qp    = qp;
        vkey.valid = true;

        const bool valueChanged = !qf || !(m_cache.vkey == vkey);
        if (valueChanged)
        {
          m_cache.vkey = vkey;

          const auto& rq = qf ? qf->getPoint(qp) : p.getReferenceCoordinates();

          // Cheap scalar P1 element (no allocations)
          const P1Element<ScalarType> feScalar(geom);
          const size_t nv = feScalar.getCount();

          for (size_t a = 0; a < nv; ++a)
            m_cache.phi_vertex[a] = feScalar.getBasis(a)(rq);

          if constexpr (std::is_same_v<RangeType, ScalarType>)
          {
            for (size_t a = 0; a < nv; ++a)
              m_cache.basis[a] = m_cache.phi_vertex[a];
          }
          else
          {
            static_assert(FormLanguage::IsVectorRange<RangeType>::Value);

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

      /**
       * @brief Gets the basis function of a local degree of freedom.
       * @param local Index in the local numbering.
       * @returns Value of the selected local basis function at the evaluation point.
       */
      constexpr
      const RangeType& getBasis(size_t local) const
      {
        assert(m_cache.skey);
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
       * @param poly Mesh entity used by this operation.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
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

namespace Rodin::Variational
{
  /// @brief Matrix-valued shape-function specialization with cached scalar basis tabulation.
  template <class Derived, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, P1<Math::SpatialMatrix<Scalar>, Mesh>, Space>
    : public ShapeFunctionBase<
        ShapeFunction<Derived, P1<Math::SpatialMatrix<Scalar>, Mesh>, Space>,
        P1<Math::SpatialMatrix<Scalar>, Mesh>, Space>
  {
    public:
      /// @brief CRTP or finite element base class.
      using FES = P1<Math::SpatialMatrix<Scalar>, Mesh>;
      /// @brief Current shape-function specialization.
      using Shape = ShapeFunction;
      /// @brief Existing shape-function interface.
      using Parent = ShapeFunctionBase<
        ShapeFunction<Derived, P1<Math::SpatialMatrix<Scalar>, Mesh>, Space>, FES, Space>;
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
