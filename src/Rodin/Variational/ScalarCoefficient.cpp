/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   void
   ScalarCoefficient<std::function<double(const double*)>>::buildMFEMCoefficient()
   {
      m_mfemCoefficient.emplace(
            [this](const mfem::Vector& v)
            {
               return m_f(v.GetData());
            });
   }

   void ScalarCoefficient<std::map<int, double>>::buildMFEMCoefficient()
   {
      int maxAttr = m_pieces.rbegin()->first;
      m_mfemCoefficient.emplace(maxAttr);
      for (int i = 1; i <= maxAttr; i++)
      {
         auto v = m_pieces.find(i);
         if (v != m_pieces.end())
            (*m_mfemCoefficient)(i) = v->second;
         else
            (*m_mfemCoefficient)(i) = 0.0;
      }
   }
}
