#ifndef RODIN_VARIATIONAL_TESTFUNCTION_H
#define RODIN_VARIATIONAL_TESTFUNCTION_H

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   template <class FEC>
   class TestFunction : public GridFunction<FEC>
   {
      private:
         TestFunction(FiniteElementSpace<FEC>& fes)
            : GridFunction<FEC>(fes)
         {}
   };
}
#endif
