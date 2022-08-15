#include "BasisOperator.h"

namespace Rodin::Variational
{
   // ---- SparseBasisOperator -----------------------------------------------
   // SparseBasisOperator(
   //       std::unique_ptr<double[]> data,
   //       const SparsityPattern& pattern,
   //       int rows, int cols, int dofs);


   // ---- BasisOperator -----------------------------------------------------
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
      for (int k = 0; k < getDOFs(); k++)
      {
         (*result)(0, 0, k) = 0.0;
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

   // mfem::DenseMatrix BasisOperator::OperatorDot(const BasisOperator& rhs) const
   // {
   //    assert(getRows() == rhs.getRows());
   //    assert(getColumns() == rhs.getColumns());
   //    mfem::DenseMatrix result(getDOFs(), rhs.getDOFs());
   //    result = 0.0;
   //    for (int i = 0; i < getDOFs(); i++)
   //    {
   //       for (int j = 0; j < rhs.getDOFs(); j++)
   //       {
   //          for (int l = 0; l < getRows(); l++)
   //          {
   //             for (int m = 0; m < getColumns(); m++)
   //             {
   //                result(i, j) += (*this)(l, m, i) * rhs(l, m, j);
   //             }
   //          }
   //       }
   //    }
   //    return result;
   // }

   void BasisOperator::addToVector(mfem::Vector& vec) const
   {
      assert(getRows() == 1 && getColumns() == 1);
      assert(getDOFs() == vec.Size());
      double* vdata = vec.GetData();
      for (int i = 0; i < getDOFs(); i++)
         vdata[i] += (*this)(0, 0, i);
   }

   // ---- SSFBO -------------------------------------------------------------
   FunctionSBO::FunctionSBO(mfem::Vector&& shape, int vdim)
      : SparseBasisOperator(vdim, 1, shape.Size() * vdim),
        m_shape(std::move(shape)),
        m_colIndex(new int (0)),
        m_rowIndex(vdim, nullptr)
   {
      assert(vdim > 0);

      for (int i = 0; i < vdim; i++)
      {
         m_rowIndex[i] = new int[vdim + 1];
         std::fill_n(m_rowIndex[i], i + 1, 0);
         std::fill_n(m_rowIndex[i] + i + 1, vdim - i, 1);
      }

      int n = m_shape.Size();
      std::vector<mfem::SparseMatrix> data;
      data.reserve(n * vdim);

      for (int i = 0; i < vdim; i++)
      {
         for (int j = 0; j < n; j++)
         {
            data.emplace_back(
                  m_rowIndex[i], m_colIndex, m_shape.GetData() + j, vdim, 1,
                  false, // Do not own index
                  false, // Do not own data
                  false);
         }
      }
      assert(n * vdim > 0);
      assert(data.size() == static_cast<size_t>(n * vdim));
      setData(std::move(data));
   }

   FunctionSBO::FunctionSBO(const FunctionSBO& other)
      : SparseBasisOperator(other),
        m_shape(other.m_shape),
        m_colIndex(new int (0)),
        m_rowIndex(other.m_rowIndex.size(), nullptr)
   {
      const int n = m_shape.Size();
      for (size_t i = 0; i < other.m_rowIndex.size(); i++)
      {
         m_rowIndex[i] = new int[n];
         std::copy(other.m_rowIndex[i], other.m_rowIndex[i] + n, m_rowIndex[i]);
      }
   }

   FunctionSBO::FunctionSBO(FunctionSBO&& other)
      : SparseBasisOperator(std::move(other)),
        m_shape(std::move(other.m_shape)),
        m_colIndex(other.m_colIndex),
        m_rowIndex(std::move(other.m_rowIndex))
   {
      other.m_colIndex = nullptr;
   }

   FunctionSBO::~FunctionSBO()
   {
      if (m_colIndex)
         delete m_colIndex;
      for (auto& v : m_rowIndex)
         delete[] v;
   }

   // ---- JSSFBO ------------------------------------------------------------
   JacobianSBO::JacobianSBO(mfem::DenseMatrix&& dshape, int sdim, int vdim)
      :  SparseBasisOperator(sdim, vdim, dshape.Height() * vdim),
         m_sdim(sdim),
         m_vdim(vdim),
         m_colIndex(vdim, nullptr),
         m_rowIndex(new int[sdim + 1])
   {
      assert(vdim > 0);

      // TODO: JacobianSBO expects the matrix to be in row-major order.
      // Hence we call transpose on the matrix. However, we should be
      // able to optimize this call away.
      dshape.Transpose();

      assert(dshape.Height() == sdim);

      // "Move" the DenseMatrix
      m_dshape.SetSize(dshape.Height(), dshape.Width());
      m_dshape.GetMemory() = std::move(dshape.GetMemory());
      dshape.Reset(nullptr, 0, 0);

      // Build the indices
      int n = m_dshape.Width();

      std::iota(m_rowIndex, m_rowIndex + sdim + 1, 0);
      for (int i = 0; i < vdim; i++)
      {
         m_colIndex[i] = new int[sdim];
         std::fill(m_colIndex[i], m_colIndex[i] + sdim, i);
      }

      std::vector<mfem::SparseMatrix> data;
      data.reserve(n * vdim);

      for (int i = 0; i < vdim; i++)
      {
         for (int j = 0; j < n; j++)
         {
            data.emplace_back(
                  m_rowIndex, m_colIndex[i], m_dshape.GetData() + sdim * j, sdim, vdim,
                  false, // Do not own index
                  false, // Do not own data
                  false);
         }
      }
      assert(n * vdim > 0);
      assert(data.size() == static_cast<size_t>(n * vdim));

      setData(std::move(data));
   }

   JacobianSBO::~JacobianSBO()
   {
      if (m_rowIndex)
         delete[] m_rowIndex;
      for (auto& v : m_colIndex)
         delete[] v;
   }
}
