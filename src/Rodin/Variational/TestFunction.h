#ifndef RODIN_VARIATIONAL_TESTFUNCTION_H
#define RODIN_VARIATIONAL_TESTFUNCTION_H

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include "ShapeFunction.h"

namespace Rodin::Variational
{
   template <class FEC, class Trait>
   class TestFunction : public ShapeFunction<FEC, Test>
   {
      public:
         TestFunction(FiniteElementSpace<FEC, Trait>& fes)
            : ShapeFunction<FEC, Test>(fes)
         {}

         TestFunction(const TestFunction& other)
            : ShapeFunction<FEC, Test>(other)
         {}

         TestFunction(TestFunction&& other)
            : ShapeFunction<FEC, Test>(std::move(other))
         {}

         void operator=(const TestFunction&) = delete;

         void operator=(TestFunction&&) = delete;

         ShapeFunctionBase<Test>& getLeaf() override
         {
            return *this;
         }

         const ShapeFunctionBase<Test>& getLeaf() const override
         {
            return *this;

         }

         TestFunction* copy() const noexcept override
         {
            return new TestFunction(*this);
         }
   };
   template <class FEC, class Trait>
   TestFunction(FiniteElementSpace<FEC, Trait>& fes) -> TestFunction<FEC, Trait>;
}
#endif
