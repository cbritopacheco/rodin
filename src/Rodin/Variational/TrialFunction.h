#ifndef RODIN_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_VARIATIONAL_TRIALFUNCTION_H

#include "ShapeFunction.h"
#include "GridFunction.h"
#include "FiniteElementSpace.h"

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
   class TrialFunction : public GridFunction<FEC>, public TrialFunctionBase
   {
      public:
         TrialFunction(FiniteElementSpace<FEC>& fes)
            :  GridFunction<FEC>(fes)
         {}

         TrialFunction(const TrialFunction& other) = default;

         void getValue(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               ScalarShape& values) const override
         {
            fe.CalcPhysShape(trans, values);
         }

         FiniteElementSpace<FEC>& getFiniteElementSpace() override
         {
            return GridFunction<FEC>::getFiniteElementSpace();
         }

         const FiniteElementSpace<FEC>& getFiniteElementSpace() const override
         {
            return GridFunction<FEC>::getFiniteElementSpace();
         }

         ShapeFunction::ValueType getValueType() const override
         {
            if (getFiniteElementSpace().getRangeDimension() == 1)
               return ShapeFunction::Scalar;
            else
               return ShapeFunction::Vector;
         }

         TrialFunction* copy() const noexcept override
         {
            return new TrialFunction(*this);
         }
   };
}
#endif
