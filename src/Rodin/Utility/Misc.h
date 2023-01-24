#ifndef RODIN_UTILITY_MISC_H
#define RODIN_UTILITY_MISC_H

#include <set>
#include <mfem.hpp>

namespace Rodin::Utility
{
  mfem::Array<int> set2marker(const std::set<int>& s, int size);
}

#endif
