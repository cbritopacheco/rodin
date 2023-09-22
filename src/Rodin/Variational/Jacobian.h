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
#include "TensorBasis.h"
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

  template <class Derived, class Operand>
  class JacobianBase;

  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobianient of a P1 GridFunction
   */
  template <class Derived, class FESType>
  class JacobianBase<Derived, GridFunction<FESType>>
    : public MatrixFunctionBase<JacobianBase<Derived, GridFunction<FESType>>>
  {
    public:
      using FES = FESType;

      using Operand = GridFunction<FES>;

      /// Parent class
      using Parent = MatrixFunctionBase<JacobianBase<Derived, Operand>>;

      /**
       * @brief Constructs the Jacobianient of an @f$ \mathbb{P}^1 @f$ function
       * @f$ u @f$.
       * @param[in] u P1 GridFunction
       */
      JacobianBase(const Operand& u)
        : m_u(u)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      }

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
      Math::Matrix getValue(const Geometry::Point& p) const
      {
        Math::Matrix out;
        const auto& polytope = p.getPolytope();
        const auto& polytopeMesh = polytope.getMesh();
        const auto& gf = m_u.get();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();
        if (polytopeMesh.isSubMesh())
        {
          const auto& submesh = polytopeMesh.asSubMesh();
          assert(submesh.getParent() == fes.getMesh());
          interpolate(out, submesh.inclusion(p));
        }
        else if (fesMesh.isSubMesh())
        {
          const auto& submesh = fesMesh.asSubMesh();
          assert(submesh.getParent() == polytopeMesh);
          interpolate(out, submesh.restriction(p));
        }
        else
        {
          interpolate(out, p);
        }
        return out;
      }

      inline
      void getValue(Math::Matrix& out, const Geometry::Point& p) const
      {
        interpolate(out, p);
      }

      inline
      constexpr
      const Operand& getOperand() const
      {
        return m_u.get();
      }

      inline
      constexpr
      auto interpolate(Math::Matrix& out, const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).interpolate(out, p);
      }

      inline
      JacobianBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const Operand> m_u;
  };
}

#endif
