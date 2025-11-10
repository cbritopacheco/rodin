/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Div.h
 * @brief Divergence operator for vector-valued functions.
 *
 * This file defines the Div class, which computes the divergence of
 * vector-valued functions in variational formulations. The divergence is
 * a fundamental differential operator that measures the "outflow" of a vector field.
 *
 * ## Mathematical Foundation
 * For a vector field @f$ \mathbf{u} : \Omega \subset \mathbb{R}^d \to \mathbb{R}^d @f$,
 * the divergence is defined as:
 * @f[
 *   \nabla \cdot \mathbf{u} = \sum_{i=1}^d \frac{\partial u_i}{\partial x_i}
 * @f]
 *
 * ## Applications
 * - Incompressibility constraint: @f$ \nabla \cdot \mathbf{u} = 0 @f$
 * - Conservation laws: @f$ \nabla \cdot \mathbf{F} = 0 @f$
 * - Mixed formulations for elliptic problems
 */
#ifndef RODIN_VARIATIONAL_DIV_H
#define RODIN_VARIATIONAL_DIV_H

#include "ForwardDecls.h"

#include "Jacobian.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
    * @defgroup DivSpecializations Div Template Specializations
    * @brief Template specializations of the Div class.
    * @see Div
    */

  /**
   * @ingroup RodinVariational
   * @brief Base class for divergence operator implementations.
   *
   * DivBase provides the foundation for computing divergences of vector-valued
   * functions. The divergence operator maps vector fields to scalar fields.
   *
   * @tparam Operand Type of the vector function
   * @tparam Derived Derived class (CRTP pattern)
   */
  template <class Operand, class Derived>
  class DivBase;

  /**
   * @ingroup DivSpecializations
   * @brief Divergence of a P1 GridFunction
   */
  template <class FES, class Data, class Derived>
  class DivBase<GridFunction<FES, Data>, Derived>
    : public ScalarFunctionBase<typename FormLanguage::Traits<FES>::ScalarType, DivBase<GridFunction<FES, Data>, Derived>>
  {
    public:
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using OperandType = GridFunction<FES, Data>;

      /// Parent class
      using Parent = ScalarFunctionBase<ScalarType, DivBase<OperandType, Derived>>;

      /**
       * @brief Constructs the Div of a @f$ \mathbb{P}_1 @f$ function @f$ u
       * @f$.
       * @param[in] u P1 GridFunction
       */
      DivBase(const OperandType& u)
        : m_u(u)
      {}

      /**
       * @brief Copy constructor
       */
      DivBase(const DivBase& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor
       */
      DivBase(DivBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      /**
       * @brief Evaluates the divergence at a point.
       * @param[in] p Point at which to evaluate
       * @return Divergence value @f$ \nabla \cdot \mathbf{u}(p) @f$
       *
       * Computes the divergence by summing the partial derivatives:
       * @f$ \nabla \cdot \mathbf{u} = \sum_{i=1}^d \frac{\partial u_i}{\partial x_i} @f$
       * Handles mesh inclusion and submesh restrictions automatically.
       */
      ScalarType getValue(const Geometry::Point& p) const
      {
        static thread_local ScalarType s_out;
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
       * @brief Gets the operand grid function.
       * @return Reference to the vector-valued grid function
       */
      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Interpolates the divergence at a point (to be overridden in derived class).
       * @param[out] out Output scalar for divergence result
       * @param[in] p Point at which to interpolate
       *
       * This virtual function is overridden in derived classes (e.g., P1::Div)
       * to provide finite element-specific divergence computation.
       */
      constexpr
      void interpolate(ScalarType& out, const Geometry::Point& p) const
      {
        static_cast<const Derived&>(*this).interpolate(out, p);
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       */
      DivBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
  };
}

#endif
