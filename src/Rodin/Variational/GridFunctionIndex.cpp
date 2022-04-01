/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include "Rodin/Core.h"
#include "GridFunction.h"
#include "FiniteElementSpace.h"

#include "GridFunctionIndex.h"

namespace Rodin::Variational
{
   GridFunctionIndex<bool>::GridFunctionIndex(const std::deque<bool>& idx)
      : m_idx(idx)
   {}

   std::vector<int> GridFunctionIndex<bool>::getIndices() const
   {
      std::vector<int> res;
      res.reserve(m_idx.size());
      for (size_t i = 0; i < m_idx.size(); i++)
      {
         if (m_idx[i])
            res.push_back(i);
      }
      return res;
   }

   GridFunctionIndex<bool>
   GridFunctionIndex<bool>::operator||(const GridFunctionIndex<bool>& other)
   {
      assert(m_idx.size() == other.m_idx.size());
      std::deque<bool> idx;
      std::transform(
            m_idx.begin(), m_idx.end(),
            other.m_idx.begin(), std::back_inserter(idx),
            std::logical_or<bool>());
      return idx;
   }
}

namespace Rodin
{
   Variational::GridFunctionIndex<bool>
   isNaN(const Variational::GridFunction<Variational::H1>& v)
   {
      std::deque<bool> idx;
      int size = v.getHandle().Size();
      for (int i = 0; i < size; i++)
         idx.push_back(Rodin::isNaN(v.getHandle()[i]));
      return idx;
   }

   Variational::GridFunctionIndex<bool>
   isInf(const Variational::GridFunction<Variational::H1>& v)
   {
      std::deque<bool> idx;
      int size = v.getHandle().Size();
      for (int i = 0; i < size; i++)
         idx.push_back(Rodin::isInf(v.getHandle()[i]));
      return idx;
   }
}
