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
  template <class Range, class Data, class Mesh>
  struct Traits<
    Variational::Jacobian<
      Variational::GridFunction<
        Variational::P0g<Range, Mesh>, Data>>>
  {
    using FESType = Variational::P0g<Range, Mesh>;
    using OperandType = Variational::GridFunction<FESType, Data>;
  };

  template <class NestedDerived, class Range, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Jacobian<
      Variational::ShapeFunction<NestedDerived, Variational::P0g<Range, Mesh>, Space>>>
  {
    using FESType = Variational::P0g<Range, Mesh>;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
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
  class Jacobian<GridFunction<P0g<Math::Vector<Scalar>, Mesh>, Data>> final
    : public JacobianBase<
        GridFunction<P0g<Math::Vector<Scalar>, Mesh>, Data>,
        Jacobian<GridFunction<P0g<Math::Vector<Scalar>, Mesh>, Data>>>
  {
    public:
      using FESType = P0g<Math::Vector<Scalar>, Mesh>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;
      using Parent = JacobianBase<OperandType, Jacobian<OperandType>>;

      Jacobian(const OperandType& u)
        : Parent(u)
      {}

      Jacobian(const Jacobian& other)
        : Parent(other)
      {}

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
  class Jacobian<ShapeFunction<ShapeFunctionDerived, P0g<Math::Vector<Scalar>, Mesh>, Space>> final
    : public ShapeFunctionBase<
        Jacobian<ShapeFunction<ShapeFunctionDerived, P0g<Math::Vector<Scalar>, Mesh>, Space>>,
        P0g<Math::Vector<Scalar>, Mesh>,
        Space>
  {
    public:
      using FESType = P0g<Math::Vector<Scalar>, Mesh>;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      using OperandType = ShapeFunction<ShapeFunctionDerived, FESType, SpaceType>;

      using Parent =
        ShapeFunctionBase<
          Jacobian<OperandType>,
          FESType,
          SpaceType>;

      explicit Jacobian(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr)
      {}

      Jacobian(const Jacobian& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr)
      {}

      Jacobian(Jacobian&& other)
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
  Jacobian(const ShapeFunction<ShapeFunctionDerived, P0g<Math::Vector<Scalar>, Mesh>, Space>&)
    -> Jacobian<ShapeFunction<ShapeFunctionDerived, P0g<Math::Vector<Scalar>, Mesh>, Space>>;
}

#endif
