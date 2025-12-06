#ifndef RODDIN_VARIATIONAL_RE_H
#define RODDIN_VARIATIONAL_RE_H

/**
 * @file Re.h
 * @brief Real part extraction from complex-valued functions.
 */

#include "ForwardDecls.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Extracts the real part of a complex-valued function.
   *
   * For a complex function @f$ f: \Omega \to \mathbb{C} @f$ with
   * @f$ f(x) = u(x) + iv(x) @f$, this operator extracts the real part:
   * @f[
   *    \text{Re}(f)(x) = u(x) \in \mathbb{R}
   * @f]
   *
   * @tparam Operand Type of the operand function
   *
   * @see Im, ComplexFunction, Conjugate
   */
  template <class Operand>
  class Re;

  /**
   * @brief Specialization for FunctionBase operands.
   * @tparam NestedDerived Nested derived type
   */
  template <class NestedDerived>
  class Re<FunctionBase<NestedDerived>> : public RealFunctionBase<Re<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Re<FunctionBase<NestedDerived>>>;

      /**
       * @brief Constructs Re from a complex function.
       * @param[in] f Complex-valued function operand
       */
      Re(const OperandType& f)
        : m_operand(f.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Re object to copy
       */
      Re(const Re& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Re object to move from
       */
      Re(Re&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      /**
       * @brief Gets the operand function.
       * @returns Reference to the complex operand
       */
      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      /**
       * @brief Evaluates the real part at a point.
       * @param[in] p Point at which to evaluate
       * @returns Real part @f$ u(p) @f$ where @f$ f(p) = u(p) + iv(p) @f$
       */
      constexpr
      Real getValue(const Geometry::Point& p) const
      {
        return getOperand().getValue(p).real();
      }

      /**
       * @brief Creates a polymorphic copy.
       * @returns Pointer to a new Re object
       */
      Re* copy() const noexcept override
      {
        return new Re(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  /**
   * @brief CTAD for Re.
   */
  template <class NestedDerived>
  Re(const FunctionBase<NestedDerived>&) -> Re<FunctionBase<NestedDerived>>;
}

#endif
