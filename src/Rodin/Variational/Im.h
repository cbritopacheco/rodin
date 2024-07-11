#ifndef RODDIN_VARIATIONAL_IM_H
#define RODDIN_VARIATIONAL_IM_H

#include "ForwardDecls.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  template <class Operand>
  class Im;

  template <class NestedDerived>
  class Im<FunctionBase<NestedDerived>> : public RealFunctionBase<Im<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Im<FunctionBase<NestedDerived>>>;

      Im(const OperandType& f)
        : m_operand(f.copy())
      {}

      Im(const Im& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      Im(Im&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      inline
      constexpr
      Real getValue(const Geometry::Point& p) const
      {
        return getOperand().getValue(p).imag();
      }

      inline Im* copy() const noexcept override
      {
        return new Im(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived>
  Im(const FunctionBase<NestedDerived>&) -> Im<FunctionBase<NestedDerived>>;

}

#endif
