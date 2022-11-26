#include "Assembly.h"

namespace Rodin::Variational
{
   std::ostream& operator<<(std::ostream& os, const Backend& backend)
   {
      switch (backend)
      {
         case Backend::Native:
         {
            os << "Native";
            break;
         }
      }
      return os;
   }
}
