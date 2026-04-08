#ifndef RODIN_VARIATIONAL_IM_H
#define RODIN_VARIATIONAL_IM_H

/**
 * @file Im.h
 * @brief Imaginary part extraction from complex-valued functions.
 */

#include "ForwardDecls.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Extracts the imaginary part of a complex-valued function.
   *
   * For a complex function @f$ f: \Omega \to \mathbb{C} @f$ with
   * @f$ f(x) = u(x) + iv(x) @f$, this operator extracts the imaginary part:
   * @f[
   *    \text{Im}(f)(x) = v(x) \in \mathbb{R}
   * @f]
   *
   * @tparam Operand Type of the operand function
   *
   * @see Re, ComplexFunction, Conjugate
   */
  template <class Operand>
  class Im;

  /**
   * @brief Specialization for FunctionBase operands.
   * @tparam NestedDerived Nested derived type
   */
  template <class NestedDerived>
  class Im<FunctionBase<NestedDerived>> : public RealFunctionBase<Im<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Im<FunctionBase<NestedDerived>>>;

      /**
       * @brief Constructs Im from a complex function.
       * @param[in] f Complex-valued function operand
       */
      Im(const OperandType& f)
        : m_operand(f.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Im object to copy
       */
      Im(const Im& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Im object to move from
       */
      Im(Im&& other)
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
       * @brief Evaluates the imaginary part at a point.
       * @param[in] p Point at which to evaluate
       * @returns Imaginary part @f$ v(p) @f$ where @f$ f(p) = u(p) + iv(p) @f$
       */
      constexpr
      Real getValue(const Geometry::Point& p) const
      {
        return getOperand().getValue(p).imag();
      }

      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        return GetOrderIfConstant(getOperand(), polytope);
      }

      /**
       * @brief Creates a polymorphic copy.
       * @returns Pointer to a new Im object
       */
      Im* copy() const noexcept override
      {
        return new Im(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  /**
   * @brief CTAD for Im.
   */
  template <class NestedDerived>
  Im(const FunctionBase<NestedDerived>&) -> Im<FunctionBase<NestedDerived>>;

}

#endif
