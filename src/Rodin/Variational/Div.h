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
   * @brief Base class for Div classes.
   */
  template <class Operand, class Derived>
  class DivBase;

  /**
   * @ingroup DivSpecializations
   * @brief Divergence of a P1 GridFunction
   */
  template <class FES, class Derived>
  class DivBase<GridFunction<FES>, Derived>
    : public ScalarFunctionBase<typename FormLanguage::Traits<FES>::ScalarType, DivBase<GridFunction<FES>, Derived>>
  {
    public:
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using OperandType = GridFunction<FES>;

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

      ScalarType getValue(const Geometry::Point& p) const
      {
        ScalarType out;
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
        return out;
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
