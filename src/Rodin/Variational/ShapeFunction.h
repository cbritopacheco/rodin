#ifndef RODIN_VARIATIONAL_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_SHAPEFUNCTION_H

#include <mfem.hpp>

#include "Rodin/FormLanguage/Base.h"

#include "H1.h"
#include "ForwardDecls.h"
#include "BasisOperator.h"
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

   class ShapeComputator
   {
      public:
         ShapeComputator(
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip)
            : m_trans(trans),
              m_ip(ip)
         {}

         const mfem::Vector& getPhysicalShape(const mfem::FiniteElement& el)
         {
            auto it = m_physShapeLookup.find(&el);
            if (it != m_physShapeLookup.end())
            {
               return it->second;
            }
            else
            {
               auto itt =  m_physShapeLookup.insert(it, {&el, mfem::Vector(el.GetDof())});
               el.CalcPhysShape(m_trans, itt->second);
               return itt->second;
            }
         }

         const mfem::DenseMatrix& getPhysicalDShape(const mfem::FiniteElement& el)
         {
            auto it = m_physDShapeLookup.find(&el);
            if (it != m_physDShapeLookup.end())
            {
               return it->second;
            }
            else
            {
               auto itt =  m_physDShapeLookup.insert(
                     it, {&el, mfem::DenseMatrix(el.GetDof(), m_trans.GetSpaceDim())});
               el.CalcPhysDShape(m_trans, itt->second);
               return itt->second;
            }
         }

         mfem::ElementTransformation& getElementTransformation()
         {
            return m_trans;
         }

         const mfem::ElementTransformation& getElementTransformation() const
         {
            return m_trans;
         }

         const mfem::IntegrationPoint& getIntegrationPoint() const
         {
            return m_ip;
         }

      private:
         mfem::ElementTransformation& m_trans;
         const mfem::IntegrationPoint& m_ip;

         std::map<const mfem::FiniteElement*, mfem::Vector> m_physShapeLookup;
         std::map<const mfem::FiniteElement*, mfem::DenseMatrix> m_physDShapeLookup;
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

         virtual void getOperator(
               DenseBasisOperator& op,
               const mfem::FiniteElement& fe,
               ShapeComputator& comp) const = 0;

         Transpose<ShapeFunctionBase<Space>> T() const
         {
            return Transpose<ShapeFunctionBase<Space>>(*this);
         }

         virtual FiniteElementSpaceBase& getFiniteElementSpace() = 0;

         virtual const FiniteElementSpaceBase& getFiniteElementSpace() const = 0;

         virtual ShapeFunctionBase<Space>* copy() const noexcept override = 0;
   };

   template <ShapeFunctionSpaceType Space, class Trait>
   class ShapeFunction<H1<Trait>, Space> : public ShapeFunctionBase<Space>
   {
      public:
         ShapeFunction(H1<Trait>& fes)
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

         H1<Trait>& getFiniteElementSpace() override
         {
            return m_fes;
         }

         const H1<Trait>& getFiniteElementSpace() const override
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

         void getOperator(
               DenseBasisOperator& op,
               const mfem::FiniteElement& fe,
               ShapeComputator& comp
               ) const override
         {
            const auto& shape = comp.getPhysicalShape(fe);
            const int n = shape.Size();
            const int vdim = getFiniteElementSpace().getVectorDimension();
            op.setSize(vdim, 1, vdim * n);
            op = 0.0;
            for (int i = 0; i < vdim; i++)
               for (int j = 0; j < n; j++)
                  op(i, 0, j + i * n) = shape(j);
         }

         virtual const ShapeFunction<H1<Trait>, Space>& getLeaf() const override = 0;

         virtual ShapeFunction* copy() const noexcept override = 0;

      private:
         H1<Trait>& m_fes;
   };
}

#endif
