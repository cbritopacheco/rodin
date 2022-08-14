#include "BasisOperator.h"

namespace Rodin::Variational
{
   // ---- SparseBasisOperator -----------------------------------------------
   // SparseBasisOperator(
   //       std::unique_ptr<double[]> data,
   //       const SparsityPattern& pattern,
   //       int rows, int cols, int dofs);


   // ---- BasisOperator -----------------------------------------------------
   std::unique_ptr<BasisOperator> BasisOperator::Trace() const
   {
      assert(getRows() == getColumns());
      auto result = new DenseBasisOperator(1, 1, getDOFs());
      (*result) = 0;
      for (int k = 0; k < getDOFs(); k++)
      {
         for (int i = 0; i < getRows(); i++)
         {
            (*result)(0, 0, k) += (*this)(i, i, k);
         }
      }
      return std::unique_ptr<BasisOperator>(result);
   }

   std::unique_ptr<BasisOperator> BasisOperator::Transpose() const
   {
      auto result = new DenseBasisOperator(getColumns(), getRows(), getDOFs());
      for (int i = 0; i < getDOFs(); i++)
      {
         for (int j = 0; j < getRows(); j++)
         {
            for (int k = 0; k < getColumns(); k++)
            {
               (*result)(k, j, i) = (*this)(j, k, i);
            }
         }
      }
      return std::unique_ptr<BasisOperator>(result);
   }

   std::unique_ptr<BasisOperator> BasisOperator::OperatorSum(const BasisOperator& rhs) const
   {
      assert(getRows() == rhs.getRows());
      assert(getColumns() == rhs.getColumns());
      assert(getDOFs() == rhs.getDOFs());
      auto result = new DenseBasisOperator(getRows(), getColumns(), getDOFs());
      for (int k = 0; k < getDOFs(); k++)
      {
         for (int i = 0; i < getRows(); i++)
         {
            for (int j = 0; j < getColumns(); j++)
            {
               (*result)(i, j, k) = (*this)(i, j, k) + rhs(i, j, k);
            }
         }
      }
      return std::unique_ptr<BasisOperator>(result);
   }

   std::unique_ptr<BasisOperator> BasisOperator::ScalarMatrixMult(
         const mfem::DenseMatrix& lhs) const
   {
      assert(getRows() == 1);
      assert(getColumns() == 1);
      auto result = new DenseBasisOperator(lhs.NumRows(), lhs.NumCols(), getDOFs());
      for (int i = 0; i < getDOFs(); i++)
      {
         (*result)(i) = lhs;
         (*result)(i) *= (*this)(0, 0, i);
      }
      return std::unique_ptr<BasisOperator>(result);
   }

   std::unique_ptr<BasisOperator> BasisOperator::LeftMatrixMult(
         const mfem::DenseMatrix& lhs) const
   {
      assert(lhs.NumCols() == getRows());
      int n = getRows();
      auto result = new DenseBasisOperator(lhs.NumRows(), getColumns(), getDOFs());
      for (int l = 0; l < getDOFs(); l++)
      {
         for (int i = 0; i < lhs.NumRows(); i++)
         {
            for (int j = 0; j < getColumns(); j++)
            {
               for (int k = 0; k < n; k++)
               {
                  (*result)(i, j, l) = lhs(i, k) * (*this)(k, j, l);
               }
            }
         }
      }
      return std::unique_ptr<BasisOperator>(result);
   }

   std::unique_ptr<BasisOperator> BasisOperator::RightMatrixMult(
         const mfem::DenseMatrix& rhs) const
   {
      assert(getColumns() == rhs.NumRows());
      int n = getColumns();
      auto result = new DenseBasisOperator(getRows(), rhs.NumCols(), getDOFs());
      for (int l = 0; l < getDOFs(); l++)
      {
         for (int i = 0; i < getRows(); i++)
         {
            for (int j = 0; j < rhs.NumCols(); j++)
            {
               for (int k = 0; k < n; k++)
               {
                  (*result)(i, j, l) = (*this)(i, k, l) * rhs(k, j);
               }
            }
         }
      }
      return std::unique_ptr<BasisOperator>(result);
   }

   std::unique_ptr<BasisOperator> BasisOperator::MatrixDot(
         const mfem::DenseMatrix& rhs) const
   {
      assert(getRows() == rhs.NumRows());
      assert(getColumns() == rhs.NumCols());
      auto result = new DenseBasisOperator(1, 1, getDOFs());
      (*result) = 0.0;
      for (int k = 0; k < getDOFs(); k++)
      {
         for (int i = 0; i < getRows(); i++)
         {
            for (int j = 0; j < getColumns(); j++)
            {
               (*result)(0, 0, k) += (*this)(i, j, k) * rhs(i, j);
            }
         }
      }
      return std::unique_ptr<BasisOperator>(result);
   }

   mfem::DenseMatrix BasisOperator::OperatorDot(const BasisOperator& rhs) const
   {
      assert(getRows() == rhs.getRows());
      assert(getColumns() == rhs.getColumns());
      mfem::DenseMatrix result(getDOFs(), rhs.getDOFs());
      result = 0.0;
      for (int i = 0; i < getDOFs(); i++)
      {
         for (int j = 0; j < rhs.getDOFs(); j++)
         {
            for (int l = 0; l < getRows(); l++)
            {
               for (int m = 0; m < getColumns(); m++)
               {
                  result(i, j) += (*this)(l, m, i) * rhs(l, m, j);
               }
            }
         }
      }
      return result;
   }

   void BasisOperator::addToVector(mfem::Vector& vec) const
   {
      assert(getRows() == 1 && getColumns() == 1);
      assert(getDOFs() == vec.Size());
      double* vdata = vec.GetData();
      for (int i = 0; i < getDOFs(); i++)
         vdata[i] += (*this)(0, 0, i);
   }

   // ---- JSSFBO ------------------------------------------------------------
   int JSSFBO::getRows() const
   {
      return m_sdim;
   }

   int JSSFBO::getColumns() const
   {
      return m_vdim;
   }

   int JSSFBO::getDOFs() const
   {
      return m_dshape.NumRows() * m_vdim;
   }

   double JSSFBO::operator()(int row, int col, int dof) const
   {
      assert(0 <= row && row < getRows());
      assert(0 <= col && col < getColumns());
      assert(0 <= dof && dof < getDOFs());
      return (static_cast<int>(dof / m_dshape.NumRows()) == col
            ) * m_dshape(dof % m_dshape.NumRows(), row);
   }

   JSSFBO& JSSFBO::operator*=(double s)
   {
      m_dshape *= s;
      return *this;
   }

   JSSFBO& JSSFBO::operator=(double s)
   {
      m_dshape = s;
      return *this;
   }

   std::unique_ptr<BasisOperator> JSSFBO::Trace() const
   {
      assert(getRows() == getColumns());
      auto result = new DenseBasisOperator(1, 1, getDOFs());
      for (int k = 0; k < getDOFs(); k++)
         (*result)(0, 0, k) = m_dshape.GetData()[k];
      return std::unique_ptr<BasisOperator>(result);
   }

   // ---- SSFBO -------------------------------------------------------------
   int SSFBO::getRows() const
   {
      return m_vdim;
   }

   int SSFBO::getColumns() const
   {
      return 1;
   }

   int SSFBO::getDOFs() const
   {
      return m_shape.Size() * m_vdim;
   }

   double SSFBO::operator()(int row, int col, int dof) const
   {
      assert(0 <= row && row < getRows());
      assert(col == 0);
      assert(0 <= dof && dof < getDOFs());
      return (static_cast<int>(dof / m_shape.Size()) == row) * m_shape(dof % m_shape.Size());
   }

   SSFBO&
   SSFBO::operator*=(double s)
   {
      m_shape *= s;
      return *this;
   }

   SSFBO&
   SSFBO::operator=(double s)
   {
      m_shape = s;
      return *this;
   }
}
