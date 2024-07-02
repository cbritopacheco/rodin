/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_JACOBIAN_H
#define RODIN_VARIATIONAL_JACOBIAN_H

#include "ForwardDecls.h"
#include "GridFunction.h"
#include "ShapeFunction.h"
#include "VectorFunction.h"
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
  template <class Derived, class Operand>
  class JacobianBase;

  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of a P1 GridFunction
   */
  template <class Derived, class FES>
  class JacobianBase<Derived, GridFunction<FES>>
    : public MatrixFunctionBase<JacobianBase<Derived, GridFunction<FES>>>
  {
    public:
      using FESType = FES;

      using OperandType = GridFunction<FES>;

      /// Parent class
      using Parent = MatrixFunctionBase<JacobianBase<Derived, OperandType>>;

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

      inline
      constexpr
      size_t getRows() const
      {
        return getOperand().getFiniteElementSpace().getVectorDimension();
      }

      inline
      constexpr
      size_t getColumns() const
      {
        return getOperand().getFiniteElementSpace().getMesh().getSpaceDimension();
      }

      inline
      Math::SpatialMatrix<Scalar> getValue(const Geometry::Point& p) const
      {
        Math::SpatialMatrix<Scalar> out;
        getValue(out, p);
        return out;
      }

      inline
      void getValue(Math::SpatialMatrix<Scalar>& out, const Geometry::Point& p) const
      {
        out.setConstant(NAN);
        const auto& polytope = p.getPolytope();
        const auto& polytopeMesh = polytope.getMesh();
        const auto& gf = getOperand();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();
        if (polytopeMesh == fesMesh)
        {
          interpolate(out, p);
        }
        else if (const auto inclusion = fesMesh.inclusion(p))
        {
          interpolate(out, *inclusion);
        }
        else if (fesMesh.isSubMesh())
        {
          const auto& submesh = fesMesh.asSubMesh();
          const auto restriction = submesh.restriction(p);
          interpolate(out, *restriction);
        }
        else
        {
          assert(false);
        }
      }

      inline
      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Interpolation function to be overriden in Derived type.
       */
      inline
      constexpr
      auto interpolate(Math::SpatialMatrix<Scalar>& out, const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).interpolate(out, p);
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       */
      inline
      JacobianBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
  };
}

#endif
