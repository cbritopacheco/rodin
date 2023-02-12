#ifndef RODIN_UTILITY_MFEM_H
#define RODIN_UTILITY_MFEM_H

#include <set>
#include <mfem.hpp>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/DenseMatrix.h"

#include "Wrap.h"

namespace Rodin::Utility
{
  inline mfem::Array<int> set2marker(const std::set<int>& s, int size)
  {
    mfem::Array<int> res(size);
    res = 0;
    for (const auto& v : s)
    {
      assert(v > 0);
      assert(v - 1 < size);
      res[v - 1] = 1;
    }
    return res;
  }
}

#endif

