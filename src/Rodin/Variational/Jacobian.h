/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_JACOBIAN_H
#define RODIN_VARIATIONAL_JACOBIAN_H

/**
 * @file
 * @brief Jacobian matrix for vector-valued functions.
 *
 * This file defines the Jacobian operator which computes the matrix of all
 * first-order partial derivatives. For a vector function
 * @f$ \mathbf{u}: \Omega \rightarrow \mathbb{R}^n @f$, the Jacobian is:
 * @f[
 *   J(\mathbf{u}) = \begin{pmatrix}
 *     \frac{\partial u_1}{\partial x_1} & \cdots & \frac{\partial u_1}{\partial x_d} \\
 *     \vdots & \ddots & \vdots \\
 *     \frac{\partial u_n}{\partial x_1} & \cdots & \frac{\partial u_n}{\partial x_d}
 *   \end{pmatrix}
 * @f]
 *
 * ## Mathematical Properties
 * - For scalar functions, the Jacobian is the gradient (row vector)
 * - The determinant @f$ \det(J) @f$ appears in change of variables
 * - The transpose @f$ J^T @f$ is often used in weak formulations
 *
 * ## Applications
 * - Nonlinear elasticity: deformation gradient
 * - Fluid mechanics: velocity gradient tensor
 * - Change of variables in integrals
 * - Newton's method for nonlinear systems
 */

#include "ForwardDecls.h"
#include "MatrixFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup JacobianSpecializations Jacobian Template Specializations
   * @brief Template specializations of the Jacobian class.
   * @see Jacobian
   */

  /**
   * @brief Base class for Jacobian classes.
   */
  template <class Operand, class Derived>
  class JacobianBase;

  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of a P1 GridFunction
   */
  template <class FES, class Data, class Derived>
  class JacobianBase<GridFunction<FES, Data>, Derived>
    : public MatrixFunctionBase<
        typename FormLanguage::Traits<FES>::ScalarType, JacobianBase<GridFunction<FES, Data>, Derived>>
  {
    public:
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;

      using Parent =
        MatrixFunctionBase<ScalarType, JacobianBase<OperandType, Derived>>;

      /**
       * @brief Constructs the Jacobianient of a @f$ \mathbb{P}_1 @f$ function
       * @f$ u @f$.
       * @param[in] u P1 GridFunction
       */
      JacobianBase(const OperandType& u)
        : m_u(u)
      {}

      /**
       * @brief Copy constructor
       */
      JacobianBase(const JacobianBase& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor
       */
      JacobianBase(JacobianBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      constexpr
      size_t getRows() const
      {
        return getOperand().getFiniteElementSpace().getVectorDimension();
      }

      constexpr
      size_t getColumns() const
      {
        return getOperand().getFiniteElementSpace().getMesh().getSpaceDimension();
      }

      decltype(auto) getValue(const Geometry::Point& p) const
      {
        static thread_local SpatialMatrixType s_out;
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

      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Interpolation function to be overriden in Derived type.
       */
      constexpr
      void interpolate(SpatialMatrixType& out, const Geometry::Point& p) const
      {
        static_cast<const Derived&>(*this).interpolate(out, p);
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       */
      JacobianBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
  };
}

#endif
