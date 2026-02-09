/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0G_GRAD_H
#define RODIN_VARIATIONAL_P0G_GRAD_H

/**
 * @file
 * @brief Gradient operator specialization for P0g (global constant) functions.
 *
 * For P0g functions (globally constant), the gradient is identically zero:
 *   ∇u = 0.
 *
 * This holds for all element types and at all points (cells/faces/boundary).
 */

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Variational/Grad.h"
#include "Rodin/Variational/ShapeFunction.h"

#include "Rodin/Variational/P0g/ForwardDecls.h"

namespace Rodin::Variational
{
  template <class Operand, class Derived>
  class GradBase;

  // ---------------------------------------------------------------------------
  // Grad of a GridFunction in P0g
  // ---------------------------------------------------------------------------
  template <class Scalar, class Mesh, class Data>
  class Grad<GridFunction<P0g<Scalar, Mesh>, Data>> final
    : public GradBase<GridFunction<P0g<Scalar, Mesh>, Data>, Grad<GridFunction<P0g<Scalar, Mesh>, Data>>>
  {
    public:
      using FESType = P0g<Scalar, Mesh>;
      using OperandType = GridFunction<FESType, Data>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using Parent = GradBase<OperandType, Grad<OperandType>>;

      explicit Grad(const OperandType& u)
        : Parent(u)
      {}

      Grad(const Grad& other)
        : Parent(other)
      {}

      Grad(Grad&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Interpolates ∇u at point p (always zero for P0g).
       */
      void interpolate(SpatialVectorType& out, const Geometry::Point& p) const
      {
        const auto& poly = p.getPolytope();
        const size_t d = poly.getDimension(); // works for cell or face
        out.resize(d);
        out.setZero();
      }

      /**
       * @brief Polynomial order of the gradient (zero).
       *
       * P0g has order 0; its gradient is identically 0.
       * Returning 0 is consistent with "zero polynomial".
       */
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }
  };

  // ---------------------------------------------------------------------------
  // Grad of a ShapeFunction in P0g
  // ---------------------------------------------------------------------------
  template <class NestedDerived, class Scalar, class Mesh, ShapeFunctionSpaceType SpaceType>
  class Grad<ShapeFunction<NestedDerived, P0g<Scalar, Mesh>, SpaceType>> final
    : public ShapeFunctionBase<Grad<ShapeFunction<NestedDerived, P0g<Scalar, Mesh>, SpaceType>>, P0g<Scalar, Mesh>, SpaceType>
  {
    public:
      using FESType = P0g<Scalar, Mesh>;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using OperandType = ShapeFunction<NestedDerived, FESType, Space>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using Parent = ShapeFunctionBase<Grad<OperandType>, FESType, Space>;

      explicit Grad(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr)
      {}

      Grad(const Grad& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr)
      {}

      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_ip(std::exchange(other.m_ip, nullptr))
      {}

      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      constexpr
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      Grad& setIntegrationPoint(const IntegrationPoint& ip)
      {
        // Keep operand aligned (even though basis is constant).
        m_u.get().setIntegrationPoint(ip);
        m_ip = &ip;

        // Dimension comes from the polytope we are integrating on.
        const auto& poly = ip.getPoint().getPolytope();
        const size_t d = poly.getDimension();

        m_zero.resize(d);
        m_zero.setZero();

        return *this;
      }

      /**
       * @brief Number of local basis functions in the gradient object.
       *
       * Same as the operand basis count: scalar P0g has 1 basis, vector P0g has vdim bases.
       */
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      /**
       * @brief Gradient of any P0g basis function is zero.
       */
      constexpr
      const SpatialVectorType& getBasis(size_t local) const
      {
        (void) local;
        assert(m_ip);
        return m_zero;
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      const IntegrationPoint* m_ip;

      // Cached zero vector sized to the current integration polytope dimension.
      SpatialVectorType m_zero;
  };
}

#endif
