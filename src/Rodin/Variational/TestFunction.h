#ifndef RODIN_VARIATIONAL_TESTFUNCTION_H
#define RODIN_VARIATIONAL_TESTFUNCTION_H

#include "ShapeFunction.h"

namespace Rodin::Variational
{
   class TestFunctionBase : public ShapeFunction
   {
      public:
         virtual SpaceType getSpaceType() const override
         {
            return Test;
         }

         virtual TestFunctionBase* copy() const noexcept override = 0;
   };

   template <class FEC>
   class TestFunction : public TestFunctionBase
   {
      public:
         TestFunction(FiniteElementSpace<FEC>& fes)
            : m_fes(fes)
         {}

         TestFunction(const TestFunction& other)
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

         TestFunction* copy() const noexcept override
         {
            return new TrialFunction(*this);
         }

      private:
         FiniteElementSpace<FEC>& m_fes;
   };
}
#endif
