#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   VectorCoefficient::VectorCoefficient(const VectorCoefficient& other)
      :  m_dimension(other.m_dimension),
         m_mfemVectorArrayCoefficient(other.m_mfemVectorArrayCoefficient)
   {
      m_values.reserve(other.m_dimension);
      for (size_t i = 0; i < m_dimension; i++)
      {
         m_values.push_back(
               std::unique_ptr<ScalarCoefficientBase>(
                  other.m_values[i]->copy()));
      }
   }

   void VectorCoefficient::buildMFEMVectorCoefficient()
   {
      for (size_t i = 0; i < m_dimension; i++)
      {
         m_values[i]->buildMFEMCoefficient();
         m_mfemVectorArrayCoefficient.Set(i, &m_values[i]->getMFEMCoefficient(), false);
      }
   }

   mfem::VectorCoefficient& VectorCoefficient::getMFEMVectorCoefficient()
   {
      return m_mfemVectorArrayCoefficient;
   }

   size_t VectorCoefficient::getDimension() const
   {
      return m_dimension;
   }

   VectorCoefficient*
   VectorCoefficient::copy() const noexcept
   {
      return new VectorCoefficient(*this);
   }

}

