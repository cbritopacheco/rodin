/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Jacobian.h
 * @brief Jacobian matrix operator for vector-valued functions.
 *
 * This file defines the Jacobian class, which computes the Jacobian matrix
 * (matrix of all first-order partial derivatives) of vector-valued functions
 * in variational formulations.
 *
 * ## Mathematical Foundation
 * For a vector-valued function @f$ \mathbf{u} : \Omega \subset \mathbb{R}^d \to \mathbb{R}^n @f$,
 * the Jacobian matrix is defined as:
 * @f[
 *   J_{ij} = \frac{\partial u_i}{\partial x_j}
 * @f]
 * resulting in an @f$ n \times d @f$ matrix.
 *
 * ## Special Cases
 * - When @f$ n = d @f$, the determinant @f$ \det(J) @f$ appears in change of variables
 * - For @f$ d = n = 2,3 @f$, related to deformation gradient in mechanics
 * - Transpose @f$ J^T @f$ gives the gradient of components
 *
 * ## Applications
 * - Nonlinear elasticity: deformation gradient tensor
 * - Fluid dynamics: velocity gradient tensor
 * - Differential geometry: metric tensor computations
 * - Coordinate transformations
 *
 * ## Usage Example
 * ```cpp
 * // Displacement gradient for elasticity
 * P1 Vh(mesh, mesh.getSpaceDimension());
 * GridFunction<P1> u(Vh);
 * auto F = Jacobian(u);  // Deformation gradient F = ∇u
 * ```
 */
#ifndef RODIN_VARIATIONAL_JACOBIAN_H
#define RODIN_VARIATIONAL_JACOBIAN_H

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
   * @ingroup RodinVariational
   * @brief Base class for Jacobian matrix operator implementations.
   *
   * JacobianBase provides the foundation for computing Jacobian matrices of
   * vector-valued functions.
   *
   * @tparam Operand Type of the vector function
   * @tparam Derived Derived class (CRTP pattern)
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
