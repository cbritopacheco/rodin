#ifndef RODIN_VARIATIONAL_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_SHAPEFUNCTION_H

#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "H1.h"
#include "ForwardDecls.h"
#include "Rank3Operator.h"
#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
   template <ShapeFunctionSpaceType Space>
   struct DualSpaceType;

   template <>
   struct DualSpaceType<Trial>
   {
      static constexpr ShapeFunctionSpaceType Value = Test;
   };

   template <>
   struct DualSpaceType<Test>
   {
      static constexpr ShapeFunctionSpaceType Value = Trial;
   };


   template <ShapeFunctionSpaceType Space>
   class ShapeFunctionBase : public FormLanguage::Base
   {
      public:
         ShapeFunctionSpaceType getSpaceType() const
         {
            return Space;
         }

         virtual int getRows(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const = 0;

         virtual int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const = 0;

         virtual int getColumns(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const = 0;

         virtual std::unique_ptr<Internal::Rank3OperatorBase> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const = 0;

         virtual FiniteElementSpaceBase& getFiniteElementSpace() = 0;

         virtual const FiniteElementSpaceBase& getFiniteElementSpace() const = 0;

         virtual ShapeFunctionBase<Space>* copy() const noexcept override = 0;
   };

   template <ShapeFunctionSpaceType Space>
   class ShapeFunction<H1, Space> : public ShapeFunctionBase<Space>
   {
      public:
         ShapeFunction(H1& fes)
            : m_fes(fes)
         {}

         ShapeFunction(const ShapeFunction& other)
            : m_fes(other.m_fes)
         {}

         ShapeFunction(ShapeFunction&& other)
            : m_fes(other.m_fes)
         {}

         H1& getFiniteElementSpace() override
         {
            return m_fes;
         }

         const H1& getFiniteElementSpace() const override
         {
            return m_fes;
         }

         int getRows(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation&) const override
         {
            return getFiniteElementSpace().getVectorDimension();
         }

         int getColumns(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation&) const override
         {
            return 1;
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation&) const
         {
            return fe.GetDof() * getFiniteElementSpace().getVectorDimension();
         }

         std::unique_ptr<Internal::Rank3OperatorBase> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            // TODO: Performance is shitty with this representation, so create
            // a new Rank3OperatorBase derived class which has the fast
            // operations.
            int dofs = fe.GetDof();
            int vdim = getFiniteElementSpace().getVectorDimension();
            mfem::Vector phi;
            phi.SetSize(dofs);
            fe.CalcPhysShape(trans, phi);
            auto result =
               new Internal::Rank3Operator(getRows(fe, trans), getColumns(fe, trans), getDOFs(fe, trans));
            (*result) = 0.0;
            for (int i = 0; i < dofs; i++)
               for (int j = 0; j < vdim; j++)
                  (*result)(j, 0, i + j * dofs) = phi(i);
            return std::unique_ptr<Internal::Rank3OperatorBase>(result);
         }

         virtual ShapeFunction* copy() const noexcept override = 0;

      private:
         H1& m_fes;
   };
}

#endif
