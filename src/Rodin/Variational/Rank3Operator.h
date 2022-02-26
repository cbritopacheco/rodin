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
    * consisting of @f$ n @f$ degrees of freedom. Then, if
    * @f A : V_h \rightarrow V_h @f$ is a
    * linear operator of dimensions @f$ p \times q @f$, the associated operator
    * @f$ T @f$ of @f$ A(u) @f$ is a rank-3 tensor of dimensions
    * @f$ p \times n \times q@f$.
    *
    * This tensor is internally handled as a multi-dimensional array:
    * @f[
    *    T = \begin{pmatrix}
    *       T_1\\
    *       \vdots\\
    *       T_p
    *    \end{bmatrix}
    *    where each @f$ T_i @f$ is a matrix of dimensions @f$ n \times q @f$.
    * @f]
    */
   class Rank3Operator
   {
      public:
         Rank3Operator()
            : m_p(0), m_dofs(0), m_q(0)
         {}

         Rank3Operator(size_t p, size_t dofs, size_t q)
            : m_p(p), m_dofs(dofs), m_q(q)
         {
            assert(p > 0);
            assert(dofs > 0);
            assert(q > 0);
            m_tensor.reserve(p);
            for (size_t i = 0; i < p; i++)
               m_tensor.emplace_back(dofs, q);
         }

         Rank3Operator(const Rank3Operator& other)
            : m_p(other.m_p), m_dofs(other.m_dofs), m_q(other.m_q),
              m_tensor(other.m_tensor)
         {}

         Rank3Operator(Rank3Operator&& other)
            : m_p(other.m_p), m_dofs(other.m_dofs), m_q(other.m_q),
              m_tensor(std::move(other.m_tensor))
         {}

         void SetSize(size_t p, size_t dofs, size_t q)
         {
            assert(p > 0);
            assert(dofs > 0);
            assert(q > 0);
            m_tensor.resize(p);
            for (size_t i = 0; i < p; i++)
               (*this)(i).SetSize(dofs, q);
         }

         size_t GetRows() const
         {
            return m_p;
         }

         size_t GetDOFs() const
         {
            return m_dofs;
         }

         size_t GetColumns() const
         {
            return m_q;
         }

         /**
          * @brief Gets the i-th matrix of the multi-dimensional array
          */
         mfem::DenseMatrix& operator()(size_t i)
         {
            assert(i < m_p);
            return m_tensor[i];
         }

         /**
          * @brief Gets the i-th matrix of the multi-dimensional array
          */
         const mfem::DenseMatrix& operator()(size_t i) const
         {
            assert(i < m_p);
            return m_tensor[i];
         }

         /**
          * @brief Gets the (i-th, j-th, k-th) element of the multi-dimensional array
          *
          * Returns the element given by @f$ A_{(i, j, k)} = (A_i)_{(j, k)} @f$
          */
         double& operator()(size_t i, size_t j, size_t k)
         {
            assert(i < m_p);
            assert(j < m_dofs);
            assert(k < m_q);
            return m_tensor[i](j, k);
         }

         const double& operator()(size_t i, size_t j, size_t k) const
         {
            assert(i < m_p);
            assert(j < m_dofs);
            assert(k < m_q);
            return m_tensor[i](j, k);
         }

         Rank3Operator& operator*=(double s)
         {
            for (size_t i = 0; i < GetRows(); i++)
               (*this)(i) *= s;
            return *this;
         }

         Rank3Operator& operator=(double s)
         {
            for (size_t i = 0; i < GetRows(); i++)
               (*this)(i) = s;
            return *this;
         }

         Rank3Operator Transpose() const
         {
            size_t p = GetRows();
            size_t n = GetDOFs();
            size_t q = GetColumns();
            Rank3Operator result(GetColumns(), GetDOFs(), GetRows());
            for (size_t i = 0; i < p; i++)
            {
               for (size_t j = 0; j < n; j++)
               {
                  for (size_t k = 0; k < q; k++)
                  {
                     result(k, j, i) = (*this)(i, j, k);
                  }
               }
            }
            return result;
         }

         Rank3Operator VectorDot(const mfem::Vector& rhs) const
         {
            mfem::DenseMatrix tmp;
            if (GetRows() == 1)
            {
               assert(GetColumns() == rhs.Size());
               tmp.SetSize(1, GetColumns());
            }
            else
            {
               assert(GetColumns() == 1);
               assert(GetRows() == rhs.Size());
               tmp.SetSize(GetRows(), 1);
            }
            tmp.GetFromVector(0, rhs);
            return MatrixDot(std::move(tmp));
         }

         Rank3Operator MatrixDot(mfem::DenseMatrix rhs) const
         {
            assert(GetRows() == rhs.NumRows());
            assert(GetColumns() == rhs.NumCols());
            rhs.Transpose();
            return RightMatrixProduct(rhs).Trace();
         }

         Rank3Operator Trace() const
         {
            assert(GetRows() == GetColumns());
            Rank3Operator result(1, GetDOFs(), 1);
            result = 0;
            for (size_t j = 0; j < GetRows(); j++)
            {
               for (size_t i = 0; i < GetDOFs(); i++)
               {
                  result(0, i, 0) += (*this)(j, i, j);
               }
            }
            return result;
         }

         /**
          * @f[
          *    \lambda A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{p \times q} @f$,
          * @f$ \lambda \in \mathbb{R} @f$.
          */
         Rank3Operator ScalarProduct(double lhs) const
         {
            Rank3Operator result(*this);
            result *= lhs;
            return result;
         }

         /**
          * @f[
          *    \Lambda A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         Rank3Operator MatrixScalarProduct(const mfem::DenseMatrix& lhs) const
         {
            assert(GetRows() == 1);
            assert(GetColumns() == 1);

            size_t p = lhs.NumRows();
            size_t n = GetDOFs();
            size_t q = lhs.NumCols();

            Rank3Operator result(p, n, q);
            for (size_t i = 0; i < p; i++)
            {
               for (size_t j = 0; j < q; j++)
               {
                  for (size_t k = 0; k < n; k++)
                  {
                     result(i, k, j) = lhs(i, j) * (*this)(0, k, 0);
                  }
               }
            }
            return result;
         }

         /**
          * @f[
          *    A(u) \Lambda
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{s \times p} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         Rank3Operator RightMatrixProduct(const mfem::DenseMatrix& rhs) const
         {
            assert(GetColumns() == rhs.NumRows());
            Rank3Operator result(GetRows(), GetDOFs(), rhs.NumCols());
            for (size_t i = 0; i < GetRows(); i++)
               mfem::Mult((*this)(i), rhs, result(i));
            return result;
         }

         /**
          * @f[
          *    \Lambda A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{q \times s} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         Rank3Operator LeftMatrixProduct(const mfem::DenseMatrix& lhs) const
         {
            assert(GetRows() == lhs.NumCols());
            Rank3Operator result(lhs.NumRows(), GetDOFs(), GetColumns());
            for (size_t i = 0; i < GetRows(); i++)
               mfem::Mult(lhs, (*this)(i), result(i));
            return result;
         }

         /**
          * @f$ u @f$ with @f$ n @f$ degrees
          * @f$ v @f$ with @f$ v @f$ degrees
          * @f[
          *    A(u) : B(v) \rightarrow \sum V_j U_j^t
          * @f]
          * Stiffness matrix of @f$ m \times n @f$.
          */
         mfem::DenseMatrix OperatorDot(const Rank3Operator& rhs) const
         {
            mfem::DenseMatrix result(rhs.GetDOFs(), GetDOFs());
            result = 0.0;

            mfem::DenseMatrix mult(rhs.GetDOFs(), GetDOFs());
            for (size_t i = 0; i < GetRows(); i++)
            {
               mfem::MultABt(rhs(i), (*this)(i), mult);
               result += mult;
            }
            return result;
         }
      private:
         size_t m_p, m_dofs, m_q;
         std::vector<mfem::DenseMatrix> m_tensor;
   };
}

#endif
