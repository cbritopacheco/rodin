/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Jacobian.h
 * @brief Jacobian operator specialization for P0g (global constant) vector functions.
 *
 * For a globally constant vector field u in P0g, all spatial derivatives vanish:
 *   J(u) = 0.
 *
 * This holds on cells and on faces/boundary (trace choice irrelevant since it is zero).
 */
#ifndef RODIN_VARIATIONAL_P0G_JACOBIAN_H
#define RODIN_VARIATIONAL_P0G_JACOBIAN_H

#include <type_traits>
#include <utility>

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/ShapeFunction.h"

#include "Rodin/Variational/P0g/ForwardDecls.h"

namespace Rodin::FormLanguage
{
  /// @brief Type traits for @c Jacobian over a grid function: exposes the finite element
  /// space and the operand type.
  template <class Range, class Data, class Mesh>
  struct Traits<
    Variational::Jacobian<
      Variational::GridFunction<
        Variational::P0g<Range, Mesh>, Data>>>
  {
      /// @brief Finite element space type.
      using FESType = Variational::P0g<Range, Mesh>;
      /// @brief Operand type.
      using OperandType = Variational::GridFunction<FESType, Data>;
  };

  /// @brief Type traits for @c Jacobian over a shape function: exposes the finite element
  /// space, the shape function space and the operand type.
  template <class NestedDerived, class Range, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Jacobian<
      Variational::ShapeFunction<NestedDerived, Variational::P0g<Range, Mesh>, Space>>>
  {
      /// @brief Finite element space type.
      using FESType = Variational::P0g<Range, Mesh>;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
      /// @brief Operand type.
      using OperandType = Variational::ShapeFunction<NestedDerived, FESType, Space>;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of a P0g vector GridFunction (identically zero).
   */
  template <class Data, class Mesh, class Scalar>
  class Jacobian<GridFunction<P0g<Math::SpatialVector<Scalar>, Mesh>, Data>> final
    : public JacobianBase<
        GridFunction<P0g<Math::SpatialVector<Scalar>, Mesh>, Data>,
        Jacobian<GridFunction<P0g<Math::SpatialVector<Scalar>, Mesh>, Data>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = P0g<Math::SpatialVector<Scalar>, Mesh>;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      /// @brief Small spatial matrix value type.
      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      /// @brief Operand type.
      using OperandType = GridFunction<FESType, Data>;
      /// @brief Parent class type.
      using Parent = JacobianBase<OperandType, Jacobian<OperandType>>;

      /// @brief Constructs the expression from its operand.
      Jacobian(const OperandType& u)
        : Parent(u)
      {}

      /// @brief Copy constructor.
      Jacobian(const Jacobian& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      Jacobian(Jacobian&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Interpolates J(u) at point p (always zero for P0g).
       *
       * The output matrix is sized (vdim x d), where:
       * - vdim = vector dimension of the FE space (typically mesh dim)
       * - d    = dimension of the polytope we are evaluating on (cell or face)
       */
      void interpolate(SpatialMatrixType& out, const Geometry::Point& p) const
      {
        const auto& poly = p.getPolytope();
        const size_t d = poly.getDimension();

        const auto& gf  = this->getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();

        out.resize(vdim, d);
        out.setZero();
      }

      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        // Identically zero.
        return 0;
      }

      /// @brief Creates a polymorphic copy.
      Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }
  };

  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of a P0g vector ShapeFunction (identically zero).
   *
   * This is the Jacobian of the *basis functions* of P0g; since those are constant,
   * their Jacobians are zero.
   */
  template <class ShapeFunctionDerived, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  class Jacobian<ShapeFunction<ShapeFunctionDerived, P0g<Math::SpatialVector<Scalar>, Mesh>, Space>> final
    : public ShapeFunctionBase<
        Jacobian<ShapeFunction<ShapeFunctionDerived, P0g<Math::SpatialVector<Scalar>, Mesh>, Space>>,
        P0g<Math::SpatialVector<Scalar>, Mesh>,
        Space>
  {
    public:
      /// @brief Finite element space type.
      using FESType = P0g<Math::SpatialVector<Scalar>, Mesh>;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      /// @brief Small spatial matrix value type.
      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      /// @brief Operand type.
      using OperandType = ShapeFunction<ShapeFunctionDerived, FESType, SpaceType>;

      /// @brief Parent class type.
      using Parent =
        ShapeFunctionBase<
          Jacobian<OperandType>,
          FESType,
          SpaceType>;

      /// @brief Constructs the expression from its operand.
      explicit Jacobian(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr)
      {}

      /// @brief Copy constructor.
      Jacobian(const Jacobian& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr)
      {}

      /// @brief Move constructor.
      Jacobian(Jacobian&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_ip(std::exchange(other.m_ip, nullptr))
      {}

      constexpr
      /// @brief Gets the operand function.
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      constexpr
      /// @brief Gets the global DOF indices for a polytope.
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      constexpr
      /// @brief Gets the integration point the expression is evaluated at.
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /// @brief Sets the integration point the expression is evaluated at.
      Jacobian& setIntegrationPoint(const IntegrationPoint& ip)
      {
        // keep operand aligned
        m_u.get().setIntegrationPoint(ip);
        m_ip = &ip;

        const auto& poly = ip.getPoint().getPolytope();
        const size_t d = poly.getDimension();

        const auto& fes  = this->getFiniteElementSpace();
        const size_t vdim = fes.getVectorDimension();

        m_zero.resize(vdim, d);
        m_zero.setZero();

        return *this;
      }

      /**
       * @brief Returns the Jacobian of the local basis function (always zero).
       *
       * ShapeFunction Jacobian basis is a matrix (vdim x d).
       */
      constexpr
      const SpatialMatrixType& getBasis(size_t local) const
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

      Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      const IntegrationPoint* m_ip;
      SpatialMatrixType m_zero;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Jacobian of a P0g vector GridFunction
   */
  template <class Range, class Data, class Mesh>
  Jacobian(const GridFunction<P0g<Range, Mesh>, Data>&)
    -> Jacobian<GridFunction<P0g<Range, Mesh>, Data>>;

  template <class ShapeFunctionDerived, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  Jacobian(const ShapeFunction<ShapeFunctionDerived, P0g<Math::SpatialVector<Scalar>, Mesh>, Space>&)
    -> Jacobian<ShapeFunction<ShapeFunctionDerived, P0g<Math::SpatialVector<Scalar>, Mesh>, Space>>;
}

#endif
