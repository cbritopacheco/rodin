#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   template <class FEC>
   class TrialFunction : public GridFunction<FEC>
   {
      public:
         TrialFunction(FiniteElementSpace<FEC>& fes)
            : GridFunction<FEC>(fes)
         {}
   };
}
#endif
