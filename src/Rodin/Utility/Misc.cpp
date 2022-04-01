#include "Misc.h"

namespace Rodin::Utility
{
   mfem::Array<int> set2marker(const std::set<int>& s, int size)
   {
      mfem::Array<int> res(size);
      res = 0;
      for (const auto& v : s)
      {
         assert(v - 1 < size);
         res[v - 1] = 1;
      }
      return res;
   }
}
