#ifndef RODIN_VARIATIONAL_TESTFUNCTION_H
#define RODIN_VARIATIONAL_TESTFUNCTION_H

#include "Component.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  template <class FES>
  struct Traits<Variational::TestFunction<FES>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Variational::TestSpace;
  };
}

namespace Rodin::Variational
{
  template <class FES>
  class TestFunction
    : public ShapeFunction<TestFunction<FES>, FES, TestSpace>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType Space = TestSpace;

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

      constexpr
      auto x() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 1);
        return Component(*this, 0);
      }

      constexpr
      auto y() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 2);
        return Component(*this, 1);
      }

      constexpr
      auto z() const
      {
        assert(this->getFiniteElementSpace().getVectorDimension() >= 3);
        return Component(*this, 2);
      }

      constexpr
      const TestFunction& getLeaf() const
      {
        return *this;
      }

      TestFunction* copy() const noexcept override
      {
        return new TestFunction(*this);
      }
  };

  template <class FES>
  TestFunction(const FES&) -> TestFunction<FES>;
}

#endif
