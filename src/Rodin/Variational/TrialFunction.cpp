#include "Component.h"

#include "TrialFunction.h"

namespace Rodin::Variational
{
   Component<TrialFunction<H1>> TrialFunction<H1>::x() const
   {
      return Component<TrialFunction<H1>>(*this, 0);
   }

   Component<TrialFunction<H1>> TrialFunction<H1>::y() const
   {
      return Component<TrialFunction<H1>>(*this, 1);
   }

   Component<TrialFunction<H1>> TrialFunction<H1>::z() const
   {
      return Component<TrialFunction<H1>>(*this, 2);
   }
}
