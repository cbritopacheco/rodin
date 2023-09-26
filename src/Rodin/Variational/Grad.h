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
  template <class Derived, class FESType>
  class GradBase<Derived, GridFunction<FESType>>
    : public VectorFunctionBase<GradBase<Derived, GridFunction<FESType>>>
  {
    public:
      using FES = FESType;

      using Operand = GridFunction<FES>;

      /// Parent class
      using Parent = VectorFunctionBase<GradBase<Derived, Operand>>;

      /**
       * @brief Constructs the gradient of a @f$ \mathbb{P}_1 @f$ function @f$
       * u @f$.
       * @param[in] u P1 GridFunction
       */
      GradBase(const Operand& u)
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
      Math::Vector getValue(const Geometry::Point& p) const
      {
        Math::Vector out;
        const auto& polytope = p.getPolytope();
        const auto& polytopeMesh = polytope.getMesh();
        const auto& gf = m_u.get();
        const auto& fes = gf.getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();
        if (polytope.getMesh() == fes.getMesh())
        {
          interpolate(out, p);
        }
        else
        {
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
            assert(false);
            out.setConstant(NAN);
          }
        }
        return out;
      }

      inline
      void getValue(Math::Vector& out, const Geometry::Point& p) const
      {
        interpolate(out, p);
      }

      inline
      constexpr
      const Operand& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Interpolation function to be overriden in Derived type.
       */
      inline
      constexpr
      auto interpolate(Math::Vector& out, const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).interpolate(out, p);
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
      std::reference_wrapper<const Operand> m_u;
  };
}

#endif
