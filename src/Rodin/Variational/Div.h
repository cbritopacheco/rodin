#ifndef RODIN_VARIATIONAL_DIV_H
#define RODIN_VARIATIONAL_DIV_H

#include "ForwardDecls.h"

#include "Jacobian.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  /**
    * @defgroup DivSpecializations Div Template Specializations
    * @brief Template specializations of the Div class.
    * @see Div
    */

  /**
   * @brief Base class for Div classes.
   */
  template <class Derived, class Operand>
  class DivBase;

  /**
   * @ingroup DivSpecializations
   * @brief Divergence of a P1 GridFunction
   */
  template <class Derived, class FESType>
  class DivBase<Derived, GridFunction<FESType>>
    : public ScalarFunctionBase<DivBase<Derived, GridFunction<FESType>>>
  {
    public:
      using FES = FESType;

      using Operand = GridFunction<FES>;

      /// Parent class
      using Parent = ScalarFunctionBase<DivBase<Derived, Operand>>;

      /**
       * @brief Constructs the Div of a @f$ \mathbb{P}_1 @f$ function @f$ u
       * @f$.
       * @param[in] u P1 GridFunction
       */
      DivBase(const Operand& u)
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

      inline
      Scalar getValue(const Geometry::Point& p) const
      {
        Scalar out;
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
            out = NAN;
          }
        }
        return out;
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
      auto interpolate(Scalar& out, const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).interpolate(out, p);
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       */
      inline
      DivBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const Operand> m_u;
  };
}

#endif
