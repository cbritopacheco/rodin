#include "Derivative.h"

namespace Rodin::Variational
{
   Derivative Dx(GridFunction<H1>& u)
   {
      assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      return Derivative(0, 0, u);
   }

   Derivative Dy(GridFunction<H1>& u)
   {
      assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      return Derivative(1, 0, u);
   }

   Derivative Dz(GridFunction<H1>& u)
   {
      assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      return Derivative(2, 0, u);
   }
}
