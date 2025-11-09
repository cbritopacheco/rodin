/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRAD_H
#define RODIN_VARIATIONAL_GRAD_H

/**
 * @file
 * @brief Gradient operator for scalar and vector functions.
 *
 * This file defines the gradient operator @f$ \nabla @f$ which computes spatial
 * derivatives of functions. For a scalar function @f$ u: \Omega \rightarrow \mathbb{R} @f$:
 * @f[
 *   \nabla u = \left(\frac{\partial u}{\partial x_1}, \frac{\partial u}{\partial x_2}, \ldots, \frac{\partial u}{\partial x_d}\right)
 * @f]
 *
 * For a vector function @f$ \mathbf{u}: \Omega \rightarrow \mathbb{R}^n @f$, the gradient
 * is the Jacobian matrix:
 * @f[
 *   \nabla \mathbf{u} = \begin{pmatrix}
 *     \frac{\partial u_1}{\partial x_1} & \cdots & \frac{\partial u_1}{\partial x_d} \\
 *     \vdots & \ddots & \vdots \\
 *     \frac{\partial u_n}{\partial x_1} & \cdots & \frac{\partial u_n}{\partial x_d}
 *   \end{pmatrix}
 * @f]
 *
 * ## Applications
 * - Poisson equation: @f$ -\nabla \cdot (\nabla u) = f @f$
 * - Diffusion problems
 * - Gradient-based optimization
 * - Strain tensors in mechanics
 */

#include "ForwardDecls.h"

#include "VectorFunction.h"

namespace Rodin::FormLanguage
{
  template <class FES, class Data>
  struct Traits<Variational::Grad<Variational::GridFunction<FES, Data>>>
  {
    using FESType = FES;

    using OperandType = Variational::GridFunction<FESType, Data>;

    using RangeType = Math::Vector<typename FormLanguage::Traits<FESType>::ScalarType>;
  };

  template <class NestedDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Grad<Variational::ShapeFunction<NestedDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;

    using OperandType = Variational::ShapeFunction<NestedDerived, FESType, SpaceType>;

    using RangeType = Math::Vector<typename FormLanguage::Traits<FESType>::ScalarType>;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup GradSpecializations Grad Template Specializations
   * @brief Template specializations of the Grad class.
   * @see Grad
   */

  /**
   * @brief Base class for Grad classes.
   */
  template <class Operand, class Derived>
  class GradBase;

  /**
   * @ingroup GradSpecializations
   */
  template <class FES, class Data, class Derived>
  class GradBase<GridFunction<FES, Data>, Derived>
    : public VectorFunctionBase<
        typename FormLanguage::Traits<FES>::ScalarType, GradBase<GridFunction<FES, Data>, Derived>>
  {
    public:
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;

      using Parent = VectorFunctionBase<ScalarType, GradBase<OperandType, Derived>>;

      GradBase(const OperandType& u)
        : m_u(u)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      }

      /**
       * @brief Copy constructor
       */
      GradBase(const GradBase& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor
       */
      GradBase(GradBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      constexpr
      size_t getDimension() const
      {
        return m_u.get().getFiniteElementSpace().getMesh().getSpaceDimension();
      }

      decltype(auto) getValue(const Geometry::Point& p) const
      {
        static thread_local SpatialVectorType s_out;
        const auto& polytope = p.getPolytope();
        const auto& polytopeMesh = polytope.getMesh();
        const auto& gf = getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();
        if (polytopeMesh == fesMesh)
        {
          this->interpolate(s_out, p);
        }
        else if (const auto inclusion = fesMesh.inclusion(p))
        {
          this->interpolate(s_out, *inclusion);
        }
        else if (fesMesh.isSubMesh())
        {
          const auto& submesh = fesMesh.asSubMesh();
          const auto restriction = submesh.restriction(p);
          this->interpolate(s_out, *restriction);
        }
        else
        {
          assert(false);
        }
        return s_out;
      }

      /**
       * @brief Interpolation function to be overriden in Derived type.
       */
      constexpr
      void interpolate(SpatialVectorType& out, const Geometry::Point& p) const
      {
        static_cast<const Derived&>(*this).interpolate(out, p);
      }

      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       */
      GradBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Grad of a GridFunction
   */
  template <class FES, class Data>
  Grad(const GridFunction<FES, Data>&) -> Grad<GridFunction<FES, Data>>;

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Grad of a ShapeFunction
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Grad(const ShapeFunction<NestedDerived, FES, Space>&)
    -> Grad<ShapeFunction<NestedDerived, FES, Space>>;
}

#endif
