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
         {
            assert(dimension > 0);
         }

         Normal(const Normal& other)
            :  VectorFunctionBase(other),
               m_dimension(other.m_dimension)
         {}

         Normal(Normal&& other)
            :  VectorFunctionBase(std::move(other)),
               m_dimension(other.m_dimension)
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
            assert(
               // We are on a boundary element of a d-mesh in d-space
               (
                  trans.mesh->Dimension() == trans.mesh->SpaceDimension() &&
                  trans.ElementType == mfem::ElementTransformation::BDR_ELEMENT
               ) ||
               // Or we are on an element of a d-mesh in (d + 1)-space.
               (
                  trans.mesh->Dimension() == (trans.mesh->SpaceDimension() - 1) &&
                  trans.ElementType == mfem::ElementTransformation::ELEMENT
               )
            );
            value.SetSize(m_dimension);
            mfem::CalcOrtho(trans.Jacobian(), value);
            const double norm = value.Norml2();
            assert(norm > 0.0);
            switch (trans.ElementType)
            {
               case mfem::ElementTransformation::BDR_ELEMENT:
               {
                  value /= norm * (
                        1.0 - 2.0 * trans.mesh->FaceIsInterior(
                           trans.mesh->GetBdrFace(trans.ElementNo)));
                  break;
               }
               default:
               {
                  value /= norm;
               }
            }
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
