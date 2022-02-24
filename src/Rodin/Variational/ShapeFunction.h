#ifndef RODIN_VARIATIONAL_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_SHAPEFUNCTION_H

#include <mfem.hpp>

#include "FormLanguage/Base.h"

#include "H1.h"
#include "ForwardDecls.h"
#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
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

         virtual size_t getColumns(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const = 0;

         /**
          * @brief Computes the associated operator @f$ U \in \mathbb{n \times d} @f$
          *
          * Let @f$ n @f$ denote the degrees of freedom of the basis
          * representation of @f$ u = (u_1, \ldots, u_d) @f$. Furthermore, let
          * @f$ w \in \mathbb{R^n} @f$ denote the basis weights
          * Then, the associated operator @f$ U @f$ is defined as an
          * @f$ n \times d @f$ matrix such that
          * @f[
          *    u = U^t w
          * @f]
          */
         virtual void getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               mfem::DenseMatrix& op) const = 0;

         /**
          * @brief Tensor dot product
          */
         virtual void dot(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               const mfem::DenseMatrix& M,
               mfem::DenseMatrix& UM) const
         {
            mfem::DenseMatrix U(getRows(fe, trans), getColumns(fe, trans));
            getOperator(fe, trans, U);
            mfem::Mult(U, M, UM);
         }

         /**
          * @brief Scalar-Tensor multiplication
          */
         virtual void multiplication(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               double m,
               mfem::DenseMatrix& op) const
         {
            getOperator(fe, trans, op);
            op *= m;
         }

         /**
          * @brief Tensor-Tensor multiplication
          *
          */
         virtual void multiplication(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               const mfem::DenseMatrix& M,
               mfem::DenseMatrix& UMt) const
         {
            mfem::DenseMatrix U(getRows(fe, trans), getColumns(fe, trans));
            getOperator(fe, trans, U);
            mfem::MultABt(U, M, UMt);
         }

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
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation&) const override
         {
            return fe.GetDof() * getFiniteElementSpace().getVectorDimension();
         }

         size_t getColumns(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation&) const override
         {
            return getFiniteElementSpace().getVectorDimension();
         }

         void getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               mfem::DenseMatrix& op) const override
         {
            size_t dofs = fe.GetDof();
            size_t vdim = getFiniteElementSpace().getVectorDimension();
            mfem::Vector tmp(dofs * vdim);
            fe.CalcPhysShape(trans, tmp);
            op = 0.0;
            for (size_t j = 0; j < vdim; j++)
            {
               for (size_t i = 0; i < dofs; i++)
               {
                  op(i + j * dofs, j) = tmp(i + j * dofs);
               }
            }
         }

         virtual ShapeFunction* copy() const noexcept override = 0;

      private:
         H1& m_fes;
   };
}

#endif
