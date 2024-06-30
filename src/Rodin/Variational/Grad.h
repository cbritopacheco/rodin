/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRADIENT_H
#define RODIN_VARIATIONAL_GRADIENT_H

#include "ForwardDecls.h"

#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "VectorFunction.h"

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
  template <class Derived, class Operand>
  class GradBase;

  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a P1 GridFunction
   */
  template <class Derived, class FES>
  class GradBase<Derived, GridFunction<FES>>
    : public VectorFunctionBase<GradBase<Derived, GridFunction<FES>>>
  {
    public:
      using FESType = FES;

      using OperandType = GridFunction<FESType>;

      /// Parent class
      using Parent = VectorFunctionBase<GradBase<Derived, OperandType>>;

      /**
       * @brief Constructs the gradient of a @f$ \mathbb{P}_1 @f$ function @f$
       * u @f$.
       * @param[in] u P1 GridFunction
       */
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

      inline
      constexpr
      size_t getDimension() const
      {
        return m_u.get().getFiniteElementSpace().getMesh().getSpaceDimension();
      }

      inline
      Math::SpatialVector<Scalar> getValue(const Geometry::Point& p) const
      {
        Math::SpatialVector<Scalar> out;
        getValue(out, p);
        return out;
      }

      inline
      void getValue(Math::SpatialVector<Scalar>& out, const Geometry::Point& p) const
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

      /**
       * @brief Interpolation function to be overriden in Derived type.
       */
      inline
      constexpr
      auto interpolate(Math::SpatialVector<Scalar>& out, const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).interpolate(out, p);
      }

      inline
      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       */
      inline
      GradBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
  };
}

#endif
