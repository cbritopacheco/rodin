#ifndef RODIN_VARIATIONAL_NORMAL_H
#define RODIN_VARIATIONAL_NORMAL_H

#include "ForwardDecls.h"

#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   class Normal : public VectorCoefficientBase
   {
      public:
         constexpr
         Normal(int dimension)
            : m_dimension(dimension)
         {}

         constexpr
         Normal(const Normal& other)
            : m_dimension(other.m_dimension)
         {}

         size_t getDimension() const override
         {
            return m_dimension;
         }

         void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint&) const override
         {
            assert(trans.ElementType == mfem::ElementTransformation::BDR_ELEMENT);
            value.SetSize(m_dimension);
            mfem::CalcOrtho(trans.Jacobian(), value);
         }

         Normal* copy() const noexcept override
         {
            return new Normal(*this);
         }
      private:
         const int m_dimension;
   };
}

#endif