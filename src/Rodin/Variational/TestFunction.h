#ifndef RODIN_VARIATIONAL_TESTFUNCTION_H
#define RODIN_VARIATIONAL_TESTFUNCTION_H

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include "ShapeFunction.h"

namespace Rodin::Variational
{
  template <class FES>
  class TestFunction : public ShapeFunction<FES, TestSpace>
  {
    public:
      constexpr
      TestFunction(FES& fes)
        : ShapeFunction<FES, TestSpace>(fes)
      {}

      constexpr
      TestFunction(const TestFunction& other)
        : ShapeFunction<FES, TestSpace>(other)
      {}

      constexpr
      TestFunction(TestFunction&& other)
        : ShapeFunction<FES, TestSpace>(std::move(other))
      {}

      void operator=(const TestFunction&) = delete;

      void operator=(TestFunction&&) = delete;

      const TestFunction& getLeaf() const override
      {
        return *this;

      }

      TestFunction* copy() const noexcept override
      {
        return new TestFunction(*this);
      }
  };
  template <class FES>
  TestFunction(FES& fes) -> TestFunction<FES>;
}
#endif
