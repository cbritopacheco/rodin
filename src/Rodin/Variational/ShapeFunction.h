#ifndef RODIN_VARIATIONAL_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_SHAPEFUNCTION_H

#include <mfem.hpp>

#include "Rodin/FormLanguage/Base.h"

#include "H1.h"
#include "ForwardDecls.h"
#include "Rank3Operator.h"
#include "FiniteElementSpace.h"

#include "RangeShape.h"

namespace Rodin::Variational
{
   template <ShapeFunctionSpaceType Space>
   struct DualSpaceType;

   template <>
   struct DualSpaceType<TrialSpace>
   {
      static constexpr ShapeFunctionSpaceType Value = ShapeFunctionSpaceType::Test;
   };

   template <>
   struct DualSpaceType<TestSpace>
   {
      static constexpr ShapeFunctionSpaceType Value = ShapeFunctionSpaceType::Trial;
   };

   /**
    * @brief Base class for operators
    */
   template <ShapeFunctionSpaceType Space>
   class ShapeFunctionBase : public FormLanguage::Base
   {
      public:
         ShapeFunctionBase()
            : FormLanguage::Base()
         {}

         ShapeFunctionBase(const ShapeFunctionBase& other)
            : FormLanguage::Base(other)
         {}

         ShapeFunctionBase(ShapeFunctionBase&& other)
            : FormLanguage::Base(std::move(other))
         {}

         ShapeFunctionSpaceType getSpaceType() const
         {
            return Space;
         }

         RangeType getRangeType() const
         {
            if (getRows() == 1 && getColumns() == 1)
               return RangeType::Scalar;
            else if (getRows() > 1 && getColumns() == 1)
               return RangeType::Vector;
            else
               return RangeType::Matrix;
         }

         RangeShape getRangeShape() const
         {
            return {getRows(), getColumns()};
         }

         virtual const ShapeFunctionBase<Space>& getLeaf() const = 0;

         virtual int getRows() const = 0;

         virtual int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const = 0;

         virtual int getColumns() const = 0;

         virtual std::unique_ptr<Internal::Rank3Operator> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const = 0;

         Transpose<ShapeFunctionBase<Space>> T() const
         {
            return Transpose<ShapeFunctionBase<Space>>(*this);
         }

         virtual FiniteElementSpaceBase& getFiniteElementSpace() = 0;

         virtual const FiniteElementSpaceBase& getFiniteElementSpace() const = 0;

         virtual ShapeFunctionBase<Space>* copy() const noexcept override = 0;
   };

   template <ShapeFunctionSpaceType Space>
   class ShapeFunction<H1, Space> : public ShapeFunctionBase<Space>
   {
      public:
         ShapeFunction(FiniteElementSpace<H1>& fes)
            : m_fes(fes)
         {}

         ShapeFunction(const ShapeFunction& other)
            :  ShapeFunctionBase<Space>(other),
               m_fes(other.m_fes)
         {}

         ShapeFunction(ShapeFunction&& other)
            :  ShapeFunctionBase<Space>(std::move(other)),
               m_fes(other.m_fes)
         {}

         FiniteElementSpace<H1>& getFiniteElementSpace() override
         {
            return m_fes;
         }

         const FiniteElementSpace<H1>& getFiniteElementSpace() const override
         {
            return m_fes;
         }

         int getRows() const override
         {
            return getFiniteElementSpace().getVectorDimension();
         }

         int getColumns() const override
         {
            return 1;
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation&) const override
         {
            return fe.GetDof() * getFiniteElementSpace().getVectorDimension();
         }

         std::unique_ptr<Internal::Rank3Operator> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            int dofs = fe.GetDof();
            int vdim = getFiniteElementSpace().getVectorDimension();
            mfem::Vector shape;
            shape.SetSize(dofs);
            fe.CalcPhysShape(trans, shape);
            return std::unique_ptr<Internal::Rank3Operator>(
                  new Internal::ScalarShapeR3O(std::move(shape), vdim));
         }

         virtual const ShapeFunction<H1, Space>& getLeaf() const override = 0;

         virtual ShapeFunction* copy() const noexcept override = 0;

      private:
         FiniteElementSpace<H1>& m_fes;
   };
}

#endif
