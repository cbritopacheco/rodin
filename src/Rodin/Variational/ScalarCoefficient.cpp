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

   mfem::Coefficient&
   ScalarCoefficient<std::function<double(const double*)>>::getMFEMCoefficient()
   {
      assert(m_mfemCoefficient);
      return *m_mfemCoefficient;
   }
}