#ifndef RODIN_VARIATIONAL_NORMAL_H
#define RODIN_VARIATIONAL_NORMAL_H

#include "ForwardDecls.h"

#include "VectorFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Outward unit normal.
    */
   class Normal : public VectorFunctionBase
   {
      public:
         /**
          * @brief Constructs the outward unit normal.
          */
         Normal(int dimension)
            : m_dimension(dimension)
         {}

         Normal(const Normal& other)
            : m_dimension(other.m_dimension)
         {}

         int getDimension() const override
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
            const double norm = value.Norml2();
            value /= norm * (
                  1.0 - 2.0 * trans.mesh->FaceIsInterior(trans.mesh->GetBdrFace(trans.ElementNo)));
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
