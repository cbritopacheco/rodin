#ifndef RODDIN_VARIATIONAL_RE_H
#define RODDIN_VARIATIONAL_RE_H

#include "ForwardDecls.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  template <class Operand>
  class Re;

  template <class NestedDerived>
  class Re<FunctionBase<NestedDerived>> : public RealFunctionBase<Re<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = RealFunctionBase<Re<FunctionBase<NestedDerived>>>;

      Re(const OperandType& f)
        : m_operand(f.copy())
      {}

      Re(const Re& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      Re(Re&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      constexpr
      Real getValue(const Geometry::Point& p) const
      {
        return getOperand().getValue(p).real();
      }

      Re* copy() const noexcept override
      {
        return new Re(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived>
  Re(const FunctionBase<NestedDerived>&) -> Re<FunctionBase<NestedDerived>>;
}

#endif
