#ifndef RODIN_VARIATIONAL_TESTFUNCTION_H
#define RODIN_VARIATIONAL_TESTFUNCTION_H

#include "Component.h"
#include "ShapeFunction.h"

namespace Rodin::Variational
{
  template <class FESType>
  class TestFunction final
    : public ShapeFunction<TestFunction<FESType>, FESType, TestSpace>
  {
    public:
      using FES = FESType;
      using Parent = ShapeFunction<TestFunction<FESType>, FESType, TestSpace>;

      constexpr
      TestFunction(const FES& fes)
        : Parent(fes)
      {}

      constexpr
      TestFunction(const TestFunction& other)
        : Parent(other)
      {}

      constexpr
      TestFunction(TestFunction&& other)
        : Parent(std::move(other))
      {}

      void operator=(const TestFunction&) = delete;

      void operator=(TestFunction&&) = delete;

      inline
      constexpr
      auto x() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 1);
        return Component(*this, 0);
      }

      inline
      constexpr
      auto y() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 2);
        return Component(*this, 1);
      }

      inline
      constexpr
      auto z() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 3);
        return Component(*this, 2);
      }

      inline
      constexpr
      const TestFunction& getLeaf() const
      {
        return *this;
      }

      inline TestFunction* copy() const noexcept override
      {
        return new TestFunction(*this);
      }
  };

  template <class FES>
  TestFunction(const FES&) -> TestFunction<FES>;
}

#endif
