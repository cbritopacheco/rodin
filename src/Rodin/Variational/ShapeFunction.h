#ifndef RODIN_VARIATIONAL_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_SHAPEFUNCTION_H

#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "H1.h"
#include "ForwardDecls.h"
#include "Rank3Operator.h"
#include "FiniteElementSpace.h"

namespace Rodin::Variational::Internal
{
   /**
    * @brief Optimized version for functions whose original representation
    * is of the form:
    * @f[
    *    u(x) = \left(
    *       \sum^n_{i=1} w_{1, i} \phi_i(x), \ldots, \sum^n_{i=1} w_{d, i} \phi_i(x) \right)
    * @f]
    */
   class ScalarShapeR3O : public Rank3Operator
   {
      public:
         ScalarShapeR3O(mfem::Vector shape, int vdim)
            : m_shape(shape),
              m_vdim(vdim)
         {}

         int GetRows() const override;

         int GetColumns() const override;

         int GetDOFs() const override;

         ScalarShapeR3O& operator*=(double s) override;

         ScalarShapeR3O& operator=(double s) override;

         double operator()(int row, int col, int dof) const override;

         std::unique_ptr<Rank3Operator>
         VectorDot(const mfem::Vector& rhs) const override;
      private:
         mfem::Vector m_shape;
         int m_vdim;
   };
}

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

         virtual ShapeFunctionBase<Space>& getRoot() = 0;

         virtual const ShapeFunctionBase<Space>& getRoot() const = 0;

         virtual int getRows(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const = 0;

         virtual int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const = 0;

         virtual int getColumns(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const = 0;

         virtual std::unique_ptr<Rank3Operator> getOperator(
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
               const mfem::ElementTransformation&) const override
         {
            return fe.GetDof() * getFiniteElementSpace().getVectorDimension();
         }

         std::unique_ptr<Rank3Operator> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            int dofs = fe.GetDof();
            int vdim = getFiniteElementSpace().getVectorDimension();
            mfem::Vector shape;
            shape.SetSize(dofs);
            fe.CalcPhysShape(trans, shape);
            return std::unique_ptr<Rank3Operator>(
                  new Internal::ScalarShapeR3O(std::move(shape), vdim));
         }

         virtual ShapeFunction* copy() const noexcept override = 0;

      private:
         FiniteElementSpace<H1>& m_fes;
   };
}

#endif
