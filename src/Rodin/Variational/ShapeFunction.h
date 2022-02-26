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

         virtual size_t getRows(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const = 0;

         virtual size_t getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const = 0;

         virtual size_t getColumns(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const = 0;

         virtual Internal::Rank3Operator getOperator(
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

         size_t getRows(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation&) const override
         {
            return getFiniteElementSpace().getVectorDimension();
         }

         size_t getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation&) const
         {
            return fe.GetDof() * getFiniteElementSpace().getVectorDimension();
         }

         size_t getColumns(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation&) const override
         {
            return 1;
         }

         Internal::Rank3Operator getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            Internal::Rank3Operator result(getRows(fe, trans), getDOFs(fe, trans), getColumns(fe, trans));
            result = 0.0;
            size_t dofs = fe.GetDof();
            size_t vdim = getFiniteElementSpace().getVectorDimension();
            mfem::Vector tmp(dofs);
            fe.CalcPhysShape(trans, tmp);
            for (size_t i = 0; i < vdim; i++)
            {
               for (size_t j = 0; j < dofs; j++)
               {
                  result(i, j + i * vdim, 0) = tmp(j);
               }
            }
            return result;
         }

         virtual ShapeFunction* copy() const noexcept override = 0;

      private:
         H1& m_fes;
   };
}

#endif
