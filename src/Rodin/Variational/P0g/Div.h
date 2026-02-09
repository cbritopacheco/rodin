/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0G_DIV_H
#define RODIN_VARIATIONAL_P0G_DIV_H

/**
 * @file
 * @brief Divergence operator specialization for P0g (global constant) vector functions.
 *
 * For a globally constant vector field u in P0g:
 *   div(u) = 0.
 *
 * This holds on cells and on faces/boundary (trace choice irrelevant since it is zero).
 */

#include <type_traits>
#include <utility>

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/Div.h"
#include "Rodin/Variational/ShapeFunction.h"

#include "Rodin/Variational/P0g/ForwardDecls.h"

namespace Rodin::FormLanguage
{
  template <class Scalar, class Data, class Mesh>
  struct Traits<Variational::Div<Variational::GridFunction<Variational::P0g<Math::Vector<Scalar>, Mesh>, Data>>>
  {
    using FESType = Variational::P0g<Math::Vector<Scalar>, Mesh>;
    using ScalarType = Scalar;
    using OperandType = Variational::GridFunction<FESType, Data>;
  };

  template <class NestedDerived, class Scalar, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Div<
      Variational::ShapeFunction<NestedDerived, Variational::P0g<Math::Vector<Scalar>, Mesh>, Space>>>
  {
    using FESType = Variational::P0g<Math::Vector<Scalar>, Mesh>;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
    using ScalarType = Scalar;
    using OperandType =
      Variational::ShapeFunction<NestedDerived, FESType, Space>;
  };
}

namespace Rodin::Variational
{
  template <class Operand, class Derived>
  class DivBase;

  /**
   * @ingroup DivSpecializations
   * @brief Divergence of a P0g vector GridFunction (identically zero).
   */
  template <class Scalar, class Data, class Mesh>
  class Div<GridFunction<P0g<Math::Vector<Scalar>, Mesh>, Data>> final
    : public DivBase<
        GridFunction<P0g<Math::Vector<Scalar>, Mesh>, Data>,
        Div<GridFunction<P0g<Math::Vector<Scalar>, Mesh>, Data>>>
  {
    public:
      using FESType = P0g<Math::Vector<Scalar>, Mesh>;
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using OperandType = GridFunction<FESType, Data>;
      using Parent = DivBase<OperandType, Div<OperandType>>;

      explicit Div(const OperandType& u)
        : Parent(u)
      {}

      Div(const Div& other)
        : Parent(other)
      {}

      Div(Div&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Interpolates div(u) at point p (always zero for P0g).
       */
      void interpolate(ScalarType& out, const Geometry::Point&) const
      {
        out = ScalarType(0);
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      Div* copy() const noexcept override
      {
        return new Div(*this);
      }
  };

  /**
   * @ingroup DivSpecializations
   * @brief Divergence of a P0g vector ShapeFunction (identically zero).
   */
  template <class NestedDerived, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  class Div<ShapeFunction<NestedDerived, P0g<Math::Vector<Scalar>, Mesh>, Space>> final
    : public ShapeFunctionBase<
        Div<ShapeFunction<NestedDerived, P0g<Math::Vector<Scalar>, Mesh>, Space>>,
        P0g<Math::Vector<Scalar>, Mesh>,
        Space>
  {
    public:
      using FESType = P0g<Math::Vector<Scalar>, Mesh>;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using OperandType = ShapeFunction<NestedDerived, FESType, SpaceType>;

      using Parent = ShapeFunctionBase<Div<OperandType>, FESType, SpaceType>;

      explicit Div(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr),
          m_zero(ScalarType(0))
      {}

      Div(const Div& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr),
          m_zero(other.m_zero)
      {}

      Div(Div&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_zero(std::exchange(other.m_zero, ScalarType(0)))
      {}

      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      constexpr
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      Div& setIntegrationPoint(const IntegrationPoint& ip)
      {
        // keep operand aligned
        m_u.get().setIntegrationPoint(ip);
        m_ip = &ip;
        m_zero = ScalarType(0);
        return *this;
      }

      /**
       * @brief Returns div(phi_local) (always zero).
       */
      constexpr
      ScalarType getBasis(size_t local) const
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

      Div* copy() const noexcept override
      {
        return new Div(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      const IntegrationPoint* m_ip;
      ScalarType m_zero;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Div of a P0g GridFunction
   */
  template <class Scalar, class Data, class Mesh>
  Div(const GridFunction<P0g<Math::Vector<Scalar>, Mesh>, Data>&)
    -> Div<GridFunction<P0g<Math::Vector<Scalar>, Mesh>, Data>>;

  template <class NestedDerived, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  Div(const ShapeFunction<NestedDerived, P0g<Math::Vector<Scalar>, Mesh>, Space>&)
    -> Div<ShapeFunction<NestedDerived, P0g<Math::Vector<Scalar>, Mesh>, Space>>;
}

#endif
