#include "Rank3Operator.h"

namespace Rodin::Variational::Internal
{
   std::unique_ptr<Rank3Operator> Rank3Operator::Trace() const
   {
      assert(GetRows() == GetColumns());
      auto result = new DenseRank3Operator(1, 1, GetDOFs());
      (*result) = 0;
      for (int k = 0; k < GetDOFs(); k++)
      {
         for (int i = 0; i < GetRows(); i++)
         {
            (*result)(0, 0, k) += (*this)(i, i, k);
         }
      }
      return std::unique_ptr<Rank3Operator>(result);
   }

   std::unique_ptr<Rank3Operator> Rank3Operator::Transpose() const
   {
      auto result = new DenseRank3Operator(GetColumns(), GetRows(), GetDOFs());
      for (int i = 0; i < GetDOFs(); i++)
      {
         for (int j = 0; j < GetRows(); j++)
         {
            for (int k = 0; k < GetColumns(); k++)
            {
               (*result)(k, j, i) = (*this)(j, k, i);
            }
         }
      }
      return std::unique_ptr<Rank3Operator>(result);
   }

   std::unique_ptr<Rank3Operator> Rank3Operator::OperatorSum(const Rank3Operator& rhs) const
   {
      assert(GetRows() == rhs.GetRows());
      assert(GetColumns() == rhs.GetColumns());
      assert(GetDOFs() == rhs.GetDOFs());
      auto result = new DenseRank3Operator(GetRows(), GetColumns(), GetDOFs());
      for (int k = 0; k < GetDOFs(); k++)
      {
         for (int i = 0; i < GetRows(); i++)
         {
            for (int j = 0; j < GetColumns(); j++)
            {
               (*result)(i, j, k) = (*this)(i, j, k) + rhs(i, j, k);
            }
         }
      }
      return std::unique_ptr<Rank3Operator>(result);
   }

   std::unique_ptr<Rank3Operator> Rank3Operator::ScalarMatrixMult(
         const mfem::DenseMatrix& lhs) const
   {
      assert(GetRows() == 1);
      assert(GetColumns() == 1);
      auto result = new DenseRank3Operator(lhs.NumRows(), lhs.NumCols(), GetDOFs());
      for (int i = 0; i < GetDOFs(); i++)
      {
         (*result)(i) = lhs;
         (*result)(i) *= (*this)(0, 0, i);
      }
      return std::unique_ptr<Rank3Operator>(result);
   }

   std::unique_ptr<Rank3Operator> Rank3Operator::ScalarVectorMult(
         const mfem::Vector& lhs) const
   {
      assert(GetRows() == 1);
      assert(GetColumns() == 1);
      auto result = new DenseRank3Operator(lhs.Size(), 1, GetDOFs());
      for (int l = 0; l < GetDOFs(); l++)
         for (int i = 0; i < lhs.Size(); i++)
            (*result)(i, 0, l) = lhs(i) * (*this)(0, 0, l);
      return std::unique_ptr<Rank3Operator>(result);
   }

   std::unique_ptr<Rank3Operator> Rank3Operator::LeftMatrixMult(
         const mfem::DenseMatrix& lhs) const
   {
      assert(lhs.NumCols() == GetRows());
      int n = GetRows();
      auto result = new DenseRank3Operator(lhs.NumRows(), GetColumns(), GetDOFs());
      for (int l = 0; l < GetDOFs(); l++)
      {
         for (int i = 0; i < lhs.NumRows(); i++)
         {
            for (int j = 0; j < GetColumns(); j++)
            {
               for (int k = 0; k < n; k++)
               {
                  (*result)(i, j, l) = lhs(i, k) * (*this)(k, j, l);
               }
            }
         }
      }
      return std::unique_ptr<Rank3Operator>(result);
   }

   std::unique_ptr<Rank3Operator> Rank3Operator::RightMatrixMult(
         const mfem::DenseMatrix& rhs) const
   {
      assert(GetColumns() == rhs.NumRows());
      int n = GetColumns();
      auto result = new DenseRank3Operator(GetRows(), rhs.NumCols(), GetDOFs());
      for (int l = 0; l < GetDOFs(); l++)
      {
         for (int i = 0; i < GetRows(); i++)
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
      return std::unique_ptr<Rank3Operator>(result);
   }

   std::unique_ptr<Rank3Operator> Rank3Operator::VectorDot(
         const mfem::Vector& rhs) const
   {
      assert(GetRows() == 1 || GetColumns() == 1);
      auto result = new DenseRank3Operator(1, 1, GetDOFs());
      (*result) = 0.0;
      if (GetRows() == 1)
      {
         for (int i = 0; i < GetDOFs(); i++)
            for (int j = 0; j < GetColumns(); j++)
               (*result)(0, 0, i) += (*this)(0, j, i) * rhs(j);
      }
      else
      {
         assert(GetColumns() == 1);
         for (int i = 0; i < GetDOFs(); i++)
            for (int j = 0; j < GetRows(); j++)
               (*result)(0, 0, i) += (*this)(j, 0, i) * rhs(j);
      }
      return std::unique_ptr<Rank3Operator>(result);
   }

   std::unique_ptr<Rank3Operator> Rank3Operator::MatrixDot(
         const mfem::DenseMatrix& rhs) const
   {
      assert(GetRows() == rhs.NumRows());
      assert(GetColumns() == rhs.NumCols());
      auto result = new DenseRank3Operator(1, 1, GetDOFs());
      (*result) = 0.0;
      for (int k = 0; k < GetDOFs(); k++)
      {
         for (int i = 0; i < GetRows(); i++)
         {
            for (int j = 0; j < GetColumns(); j++)
            {
               (*result)(0, 0, k) += (*this)(i, j, k) * rhs(i, j);
            }
         }
      }
      return std::unique_ptr<Rank3Operator>(result);
   }

   mfem::DenseMatrix Rank3Operator::OperatorDot(const Rank3Operator& rhs) const
   {
      assert(GetRows() == rhs.GetRows());
      assert(GetColumns() == rhs.GetColumns());
      mfem::DenseMatrix result(GetDOFs(), rhs.GetDOFs());
      result = 0.0;
      for (int i = 0; i < GetDOFs(); i++)
      {
         for (int j = 0; j < rhs.GetDOFs(); j++)
         {
            for (int l = 0; l < GetRows(); l++)
            {
               for (int m = 0; m < GetColumns(); m++)
               {
                  result(i, j) += (*this)(l, m, i) * rhs(l, m, j);
               }
            }
         }
      }
      return result;
   }

   void Rank3Operator::AddToVector(mfem::Vector& vec) const
   {
      assert(GetRows() == 1 && GetColumns() == 1);
      assert(GetDOFs() == vec.Size());
      double* vdata = vec.GetData();
      for (int i = 0; i < GetDOFs(); i++)
         vdata[i] += (*this)(0, 0, i);
   }

   int JacobianShapeR3O::GetRows() const
   {
      return m_sdim;
   }

   int JacobianShapeR3O::GetColumns() const
   {
      return m_vdim;
   }

   int JacobianShapeR3O::GetDOFs() const
   {
      return m_dshape.NumRows() * m_vdim;
   }

   double JacobianShapeR3O::operator()(int row, int col, int dof) const
   {
      assert(0 <= row && row < GetRows());
      assert(0 <= col && col < GetColumns());
      assert(0 <= dof && dof < GetDOFs());
      return (static_cast<int>(dof / m_dshape.NumRows()) == col
            ) * m_dshape(dof % m_dshape.NumRows(), row);
   }

   JacobianShapeR3O& JacobianShapeR3O::operator*=(double s)
   {
      m_dshape *= s;
      return *this;
   }

   JacobianShapeR3O& JacobianShapeR3O::operator=(double s)
   {
      m_dshape = s;
      return *this;
   }

   std::unique_ptr<Rank3Operator> JacobianShapeR3O::Trace() const
   {
      assert(GetRows() == GetColumns());
      auto result = new DenseRank3Operator(1, 1, GetDOFs());
      for (int k = 0; k < GetDOFs(); k++)
         (*result)(0, 0, k) = m_dshape.GetData()[k];
      return std::unique_ptr<Rank3Operator>(result);
   }

   int ScalarShapeR3O::GetRows() const
   {
      return m_vdim;
   }

   int ScalarShapeR3O::GetColumns() const
   {
      return 1;
   }

   int ScalarShapeR3O::GetDOFs() const
   {
      return m_shape.Size() * m_vdim;
   }

   double ScalarShapeR3O::operator()(int row, int col, int dof) const
   {
      assert(0 <= row && row < GetRows());
      assert(col == 0);
      assert(0 <= dof && dof < GetDOFs());
      return (static_cast<int>(dof / m_shape.Size()) == row) * m_shape(dof % m_shape.Size());
   }

   ScalarShapeR3O&
   ScalarShapeR3O::operator*=(double s)
   {
      m_shape *= s;
      return *this;
   }

   ScalarShapeR3O&
   ScalarShapeR3O::operator=(double s)
   {
      m_shape = s;
      return *this;
   }

   std::unique_ptr<Rank3Operator>
   ScalarShapeR3O::VectorDot(const mfem::Vector& rhs) const
   {
      assert(GetRows() == rhs.Size());
      auto result = new DenseRank3Operator(1, 1, GetDOFs());
      for (int i = 0; i < m_vdim; i++)
      {
         for (int k = 0; k < m_shape.Size(); k++)
         {
            (*result)(0, 0, k + i * m_shape.Size()) = m_shape(k) * rhs(i);
         }
      }
      return std::unique_ptr<Rank3Operator>(result);
   }
}
