#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include "ShapeFunction.h"

namespace Rodin::Variational
{
   class TrialFunctionBase : public ShapeFunction
   {
      public:
         virtual SpaceType getSpaceType() const override
         {
            return Trial;
         }

         virtual TrialFunctionBase* copy() const noexcept override = 0;
   };

   template <class FEC>
   class TrialFunction : public TrialFunctionBase
   {
      public:
         TrialFunction(FiniteElementSpace<FEC>& fes)
            : m_fes(fes)
         {}

         TrialFunction(const TrialFunction& other)
            : m_fes(other.m_fes)
         {}

         const FiniteElementSpace<FEC>& getFiniteElementSpace() const override
         {
            return m_fes;
         }

         ShapeFunction::ValueType getValueType() const override
         {
            if (getFiniteElementSpace().getRangeDimension() == 1)
               return ShapeFunction::Scalar;
            else
               return ShapeFunction::Vector;
         }

         void getValue(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               ScalarShape& values) const override
         {
            fe.CalcPhysShape(trans, values);
         }

         void getValue(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               VectorShape& values) const override
         {
            fe.CalcPhysVShape(trans, values);
         }

         TrialFunction* copy() const noexcept override
         {
            return new TrialFunction(*this);
         }

      private:
         FiniteElementSpace<FEC>& m_fes;
   };
}
#endif
