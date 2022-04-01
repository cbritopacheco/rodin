#include "Rank3Operator.h"

namespace Rodin::Variational
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
      auto result = new DenseRank3Operator(GetColumns(), GetRows(), GetDOFs());
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
}
