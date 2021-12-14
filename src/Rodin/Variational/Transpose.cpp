/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Transpose.h"

namespace Rodin::Variational
{
   Transpose::Transpose(const MatrixCoefficientBase& m)
      : m_matrix(m.copy())
   {}

   Transpose::Transpose(const Transpose& other)
      : m_matrix(other.m_matrix->copy())
   {}

   int Transpose::getRows() const
   {
      return m_matrix->getColumns();
   }

   int Transpose::getColumns() const
   {
      return m_matrix->getRows();
   }

   void Transpose::buildMFEMMatrixCoefficient()
   {
      m_matrix->buildMFEMMatrixCoefficient();
      m_mfemMatrixCoefficient.emplace(m_matrix->getMFEMMatrixCoefficient());
   }

   mfem::MatrixCoefficient& Transpose::getMFEMMatrixCoefficient()
   {
      assert(m_mfemMatrixCoefficient);
      return *m_mfemMatrixCoefficient;
   }
}
