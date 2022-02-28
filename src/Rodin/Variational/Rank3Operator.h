#ifndef RODIN_VARIATIONAL_RANK3OPERATOR_H
#define RODIN_VARIATIONAL_RANK3OPERATOR_H

#include <vector>
#include <mfem.hpp>

namespace Rodin::Variational::Internal
{
   /**
    * @brief Rank-3 Tensor for discretized basis functions
    *
    * Let @f$ u \in V_h @f$ be a function which has a basis representation
    * consisting of @f$ n @f$ degrees of freedom.
    *
    * This tensor should be regarded as a multi-dimensional array:
    * @f[
    *    T = \begin{pmatrix}
    *       T_1\\
    *       \vdots\\
    *       T_n
    *    \end{bmatrix}
    *    where each @f$ T_i @f$ is a matrix of dimensions @f$ p \times q @f$.
    * @f]
    */
   class Rank3OperatorBase
   {
      public:
         virtual ~Rank3OperatorBase() = default;

         virtual void SetSize(int rows, int cols, int dofs) = 0;

         virtual int GetRows() const = 0;

         virtual int GetColumns() const = 0;

         virtual int GetDOFs() const = 0;

         /**
          * @brief Gets the i-th matrix of the multi-dimensional array
          */
         virtual mfem::DenseMatrix& operator()(int i) = 0;

         /**
          * @brief Gets the i-th matrix of the multi-dimensional array
          */
         virtual const mfem::DenseMatrix& operator()(int i) const = 0;

         /**
          * @brief Gets the (i-th, j-th, k-th) element of the multi-dimensional array
          *
          * Returns the element given by @f$ A_{(i, j, k)} = (A_i)_{(j, k)} @f$
          */
         virtual double& operator()(int i, int j, int k) = 0;

         virtual const double& operator()(int i, int j, int k) const = 0;

         virtual Rank3OperatorBase& operator*=(double s) = 0;

         virtual Rank3OperatorBase& operator=(double s) = 0;

         virtual
         std::unique_ptr<Rank3OperatorBase>
         Transpose() const = 0;

         virtual
         std::unique_ptr<Rank3OperatorBase>
         VectorDot(const mfem::Vector& rhs) const = 0;

         virtual
         std::unique_ptr<Rank3OperatorBase>
         MatrixDot(const mfem::DenseMatrix& rhs) const = 0;

         virtual
         std::unique_ptr<Rank3OperatorBase> Trace() const = 0;

         /**
          * @f[
          *    A(u) \Lambda
          * @f]
          * with @f$ A(u) \in \mathbb{R} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         virtual
         std::unique_ptr<Rank3OperatorBase>
         ScalarMatrixMult(const mfem::DenseMatrix& lhs) const = 0;

         /**
          * @f[
          *    A(u) \Lambda
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{s \times p} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         virtual
         std::unique_ptr<Rank3OperatorBase>
         RightMatrixMult(const mfem::DenseMatrix& rhs) const = 0;

         /**
          * @f[
          *    \Lambda A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{q \times s} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         virtual
         std::unique_ptr<Rank3OperatorBase>
         LeftMatrixMult(const mfem::DenseMatrix& lhs) const = 0;

         /**
          * @f$ u @f$ with @f$ n @f$ degrees
          * @f$ v @f$ with @f$ v @f$ degrees
          * @f[
          *    A(u) : B(v) \rightarrow \sum V_j U_j^t
          * @f]
          * Stiffness matrix of @f$ m \times n @f$.
          */
         virtual mfem::DenseMatrix OperatorDot(const Rank3OperatorBase& rhs) const = 0;

         virtual void AddToVector(int offset, mfem::Vector& vec) const = 0;
   };

   class Rank3Operator : public Rank3OperatorBase
   {
      public:
         Rank3Operator()
         {}

         Rank3Operator(int rows, int cols, int dofs)
            : m_tensor(rows, cols, dofs)
         {
            assert(rows > 0);
            assert(dofs > 0);
            assert(cols > 0);
         }

         Rank3Operator(const Rank3Operator& other)
            : m_tensor(other.m_tensor)
         {}

         Rank3Operator(Rank3Operator&& other)
            : m_tensor(std::move(other.m_tensor))
         {}

         Rank3Operator& operator=(double s) override
         {
            for (int i = 0; i < GetDOFs(); i++)
               (*this)(i) = s;
            return *this;
         }

         Rank3Operator& operator*=(double s) override
         {
            for (int i = 0; i < GetDOFs(); i++)
               (*this)(i) *= s;
            return *this;
         }

         void SetSize(int rows, int cols, int dofs) override
         {
            m_tensor.SetSize(rows, cols, dofs);
         }

         int GetRows() const override
         {
            return m_tensor.SizeI();
         }

         int GetColumns() const override
         {
            return m_tensor.SizeJ();
         }

         int GetDOFs() const override
         {
            return m_tensor.SizeK();
         }

         /**
          * @brief Gets the i-th matrix of the multi-dimensional array
          */
         mfem::DenseMatrix& operator()(int dof) override
         {
            return m_tensor(dof);
         }

         const mfem::DenseMatrix& operator()(int dof) const override
         {
            return m_tensor(dof);
         }

         double& operator()(int row, int col, int dof) override
         {
            return m_tensor(row, col, dof);
         }

         const double& operator()(int row, int col, int dof) const override
         {
            return m_tensor(row, col, dof);
         }

         std::unique_ptr<Rank3OperatorBase> Transpose() const override
         {
            auto result = new Rank3Operator(GetColumns(), GetRows(), GetDOFs());
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
            return std::unique_ptr<Rank3OperatorBase>(result);
         }

         std::unique_ptr<Rank3OperatorBase> VectorDot(const mfem::Vector& rhs) const override
         {
            assert(GetRows() == 1 || GetColumns() == 1);
            auto result = new Rank3Operator(1, 1, GetDOFs());
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
            return std::unique_ptr<Rank3OperatorBase>(result);
         }

         std::unique_ptr<Rank3OperatorBase>
         MatrixDot(const mfem::DenseMatrix& rhs) const override
         {
            assert(GetRows() == rhs.NumRows());
            assert(GetColumns() == rhs.NumCols());
            auto result = new Rank3Operator(1, 1, GetDOFs());
            for (int i = 0; i < GetDOFs(); i++)
               (*result)(0, 0, i) = (*this)(i) * rhs;
            return std::unique_ptr<Rank3OperatorBase>(result);
         }

         std::unique_ptr<Rank3OperatorBase> Trace() const override
         {
            assert(GetRows() == GetColumns());
            auto result = new Rank3Operator(1, 1, GetDOFs());
            for (int i = 0; i < GetDOFs(); i++)
            {
               (*result)(0, 0, i) = (*this)(i).Trace();
            }
            return std::unique_ptr<Rank3OperatorBase>(result);
         }

         std::unique_ptr<Rank3OperatorBase>
         ScalarMatrixMult(const mfem::DenseMatrix& lhs) const override
         {
            assert(GetRows() == 1);
            assert(GetColumns() == 1);
            auto result = new Rank3Operator(lhs.NumRows(), lhs.NumCols(), GetDOFs());
            for (int i = 0; i < GetDOFs(); i++)
            {
               (*result)(i) = lhs;
               (*result)(i) *= (*this)(0, 0, i);
            }
            return std::unique_ptr<Rank3OperatorBase>(result);
         }

         std::unique_ptr<Rank3OperatorBase>
         RightMatrixMult(const mfem::DenseMatrix& rhs) const override
         {
            assert(GetColumns() == rhs.NumRows());
            auto result = new Rank3Operator(GetRows(), rhs.NumCols(), GetDOFs());
            for (int i = 0; i < GetDOFs(); i++)
               mfem::Mult((*this)(i), rhs, (*result)(i));
            return std::unique_ptr<Rank3OperatorBase>(result);
         }

         std::unique_ptr<Rank3OperatorBase>
         LeftMatrixMult(const mfem::DenseMatrix& lhs) const override
         {
            assert(lhs.NumCols() == GetRows());
            auto result = new Rank3Operator(lhs.NumRows(), GetColumns(), GetDOFs());
            for (int i = 0; i < GetDOFs(); i++)
               mfem::Mult(lhs, (*this)(i), (*result)(i));
            return std::unique_ptr<Rank3OperatorBase>(result);
         }

         mfem::DenseMatrix
         OperatorDot(const Rank3OperatorBase& rhs) const override
         {
            assert(GetRows() == rhs.GetRows());
            assert(GetColumns() == rhs.GetColumns());
            mfem::DenseMatrix dot(GetDOFs(), rhs.GetDOFs());
            for (int i = 0; i < GetDOFs(); i++)
            {
               for (int j = 0; j < rhs.GetDOFs(); j++)
               {
                  dot(i, j) = (*this)(i) * rhs(j);
               }
            }
            return dot;
         }

         void AddToVector(int offset, mfem::Vector& vec) const override
         {
            const int n = m_tensor.SizeI() * m_tensor.SizeJ() * m_tensor.SizeK();
            double* vdata = vec.GetData() + offset;
            for (int i = 0; i < n; i++)
               vdata[i] += m_tensor.Data()[i];
         }
      private:
         mfem::DenseTensor m_tensor;
   };
}

#endif
