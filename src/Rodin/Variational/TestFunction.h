#ifndef RODIN_VARIATIONAL_TESTFUNCTION_H
#define RODIN_VARIATIONAL_TESTFUNCTION_H

#include "ShapeFunction.h"

namespace Rodin::Variational
{
   template <>
   class TestFunction<H1> : public ShapeFunction<H1, Test>
   {
      public:
         TestFunction(H1& fes)
            : ShapeFunction<H1, Test>(fes)
         {}

         TestFunction(const TestFunction& other)
            : ShapeFunction<H1, Test>(other)
         {}

         TestFunction(TestFunction&& other)
            : ShapeFunction<H1, Test>(std::move(other))
         {}

         TestFunction* copy() const noexcept override
         {
            return new TestFunction(*this);
         }
   };
   TestFunction(H1& fes) -> TestFunction<H1>;
}
#endif
