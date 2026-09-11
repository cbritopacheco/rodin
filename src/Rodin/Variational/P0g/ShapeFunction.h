/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ShapeFunction.h
 * @brief Shape function specializations for the global-constant (P0g)
 * finite element space.
 */
#ifndef RODIN_VARIATIONAL_P0G_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_P0G_SHAPEFUNCTION_H

#include <type_traits>
#include <utility>
#include <vector>

#include "Rodin/Variational/P0g/ForwardDecls.h"
#include "Rodin/Variational/ShapeFunction.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Math/Traits.h"

namespace Rodin::Variational
{
  /// @brief Shape function expression.
  template <class Derived, class Range, class Mesh, ShapeFunctionSpaceType Space>
  class ShapeFunction<Derived, P0g<Range, Mesh>, Space>
    : public ShapeFunctionBase<
        ShapeFunction<Derived, P0g<Range, Mesh>, Space>,
        P0g<Range, Mesh>,
        Space>
  {
    public:
      /// @brief Finite element space type.
      using FESType = P0g<Range, Mesh>;
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

      static constexpr bool IsScalarRange = std::is_same_v<RangeType, ScalarType>;
      static constexpr bool IsVectorRange = FormLanguage::IsVectorRange<RangeType>::Value;

      static_assert(IsScalarRange || IsVectorRange);

      // Cache is only instantiated/used for the vector case.
      /// @brief Cached per-cell values of a vector-valued P0g shape function.
      struct VectorCache
      {
        /// @brief Vector dimension of the finite element space.
        size_t vdim = 0;
        /// @brief Cached reference basis tabulation.
        std::vector<RangeType> basis; // basis[c] = e_c (size vdim)
      };

      ShapeFunction() = delete;

      constexpr
      explicit ShapeFunction(const FESType& fes)
        : Parent(fes),
          m_ip(nullptr)
      {}

      constexpr
      /// @brief Copy constructor.
      ShapeFunction(const ShapeFunction& other)
        : Parent(other),
          m_ip(nullptr),
          m_vcache(other.m_vcache)
      {}

      constexpr
      /// @brief Move constructor.
      ShapeFunction(ShapeFunction&& other)
        : Parent(std::move(other)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_vcache(std::move(other.m_vcache))
      {}

      constexpr
      /// @brief Gets the global DOF indices for a polytope.
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
      /// @brief Gets the integration point the expression is evaluated at.
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /// @brief Sets the integration point the expression is evaluated at.
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
      /// @brief Gets the basis function of a local degree of freedom.
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
      /// @brief Returns the polynomial order used on a mesh entity.
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      constexpr
      /// @brief Gets the operand in the shape function expression.
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
