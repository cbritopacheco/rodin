/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Transpose.h"

namespace Rodin::Variational
{
   int Transpose::getRows() const
   {
      return m_matrix->getColumns();
   }

   int Transpose::getColumns() const
   {
      return m_matrix->getRows();
   }

   void Transpose::build()
   {
      m_matrix->build();
      m_mfemMatrixCoefficient.emplace(m_matrix->get());
   }

   mfem::MatrixCoefficient& Transpose::get()
   {
      assert(m_mfemMatrixCoefficient);
      return *m_mfemMatrixCoefficient;
   }
}
