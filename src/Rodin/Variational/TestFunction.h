#ifndef RODIN_VARIATIONAL_TESTFUNCTION_H
#define RODIN_VARIATIONAL_TESTFUNCTION_H

#include "ShapeFunction.h"

namespace Rodin::Variational
{
   template <class FES>
   class TestFunction : public ShapeFunction<FES, Test>
   {
      public:
         TestFunction(FES& fes)
            : ShapeFunction<FES, Test>(fes)
         {}

         TestFunction(const TestFunction& other)
            : ShapeFunction<FES, Test>(other)
         {}

         TestFunction(TestFunction&& other)
            : ShapeFunction<FES, Test>(std::move(other))
         {}

         TestFunction* copy() const noexcept override
         {
            return new TestFunction(*this);
         }
   };
   template <class FES>
   TestFunction(FES& fes) -> TestFunction<FES>;
}
#endif
