/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TRANSPOSE_H
#define RODIN_VARIATIONAL_TRANSPOSE_H

#include "ShapeFunction.h"
#include "MatrixCoefficient.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents the transpose matrix @f$ A^T @f$ of some matrix @f$ A @f$.
    *
    * For some @f$ n \times m @f$ matrix @f$ A @f$, the transpose matrix @f$
    * A^T @f$ is an @f$ m \times n @f$ matrix defined by
    * @f[
    *    {A^T}_{ij} = A_{ji}
    * @f]
    */
   template <>
   class Transpose<MatrixCoefficientBase> : public MatrixCoefficientBase
   {
      public:
         /**
          * @brief Constructs the Transpose matrix of the given matrix.
          */
         Transpose(const MatrixCoefficientBase& m)
            : m_matrix(m.copy())
         {}

         Transpose(const Transpose& other)
            :  m_matrix(other.m_matrix->copy())
         {}

         int getRows() const override
         {
            return m_matrix->getColumns();
         }

         int getColumns() const override
         {
            return m_matrix->getRows();
         }

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            m_matrix->getValue(value, trans, ip);
            value.Transpose();
         }

         Transpose* copy() const noexcept override
         {
            return new Transpose(*this);
         }

      private:
         std::unique_ptr<MatrixCoefficientBase> m_matrix;
   };
   Transpose(const MatrixCoefficientBase&)
      -> Transpose<MatrixCoefficientBase>;

   template <ShapeFunctionSpaceType Space>
   class Transpose<ShapeFunctionBase<Space>> : public ShapeFunctionBase<Space>
   {
      public:
         Transpose(const ShapeFunctionBase<Space>& op)
            : m_op(op.copy())
         {}

         Transpose(const Transpose& other)
            : m_op(other.m_op->copy())
         {}

         Transpose(Transpose&& other)
            : m_op(std::move(other))
         {}

         int getRows(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return m_op->getColumns(fe, trans);
         }

         int getColumns(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return m_op->getRows(fe, trans);
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return m_op->getDOFs(fe, trans);
         }

         std::unique_ptr<Rank3Operator> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            return m_op->getOperator(fe, trans)->Transpose();
         }

         FiniteElementSpaceBase& getFiniteElementSpace() override
         {
            return m_op->getFiniteElementSpace();
         }

         const FiniteElementSpaceBase& getFiniteElementSpace() const override
         {
            return m_op->getFiniteElementSpace();
         }

         Transpose* copy() const noexcept override
         {
            return new Transpose(*this);
         }

      private:
         std::unique_ptr<ShapeFunctionBase<Space>> m_op;
   };
}

#endif
