#ifndef RODIN_VARIATIONAL_BASISOPERATOR_H
#define RODIN_VARIATIONAL_BASISOPERATOR_H

#include <vector>
#include <mfem.hpp>

namespace Rodin::Variational
{
   /**
    * @brief Rank-3 tensor operator for discretized functions on finite bases
    *
    * Let @f$ u \in V_h @f$ be a function which has a basis representation
    * consisting of @f$ n @f$ degrees of freedom.
    *
    * This tensor should be regarded as a multi-dimensional array:
    * @f[
    *    T =
    *    \begin{bmatrix}
    *       T_1\\
    *       \vdots\\
    *       T_n
    *    \end{bmatrix}
    * @f]
    * where each @f$ T_i @f$ is a matrix of dimensions @f$ p \times q @f$.
    */
   class BasisOperator
   {
      public:
         BasisOperator() = default;

         BasisOperator(BasisOperator&&) = default;

         BasisOperator(const BasisOperator&) = default;

         BasisOperator& operator=(BasisOperator&&) = default;

         BasisOperator& operator=(const BasisOperator&) = delete;

         /**
          * @f[
          *    A(u)^T
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{p \times q} @f$.
          */
         virtual
         std::unique_ptr<BasisOperator> Transpose() const;

         /**
          * @f[
          *    \Lambda A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         virtual
         std::unique_ptr<BasisOperator> ScalarMatrixMult(const mfem::DenseMatrix& lhs) const;

         /**
          * @f[
          *    A(u) \Lambda
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{s \times p} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         virtual
         std::unique_ptr<BasisOperator> RightMatrixMult(const mfem::DenseMatrix& rhs) const;

         /**
          * @f[
          *    \Lambda A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{q \times s} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         virtual
         std::unique_ptr<BasisOperator> LeftMatrixMult(const mfem::DenseMatrix& lhs) const;

         /**
          * @f[
          *    \Lambda : A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{p \times q} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         virtual
         std::unique_ptr<BasisOperator> MatrixDot(const mfem::DenseMatrix& rhs) const;

         virtual
         std::unique_ptr<BasisOperator> OperatorSum(const BasisOperator& rhs) const;

         // virtual
         // mfem::DenseMatrix OperatorDot(const BasisOperator& rhs) const;

         virtual ~BasisOperator() = default;

         /**
          * @brief Adds to element vector.
          */
         virtual void addToVector(mfem::Vector& vec) const;

         virtual bool isDense() const = 0;

         virtual bool isSparse() const = 0;

         virtual int getRows() const = 0;

         virtual int getColumns() const = 0;

         virtual int getDOFs() const = 0;

         virtual BasisOperator& operator*=(double s) = 0;

         virtual double operator()(int row, int col, int dof) const = 0;
   };

   class DenseBasisOperator : public BasisOperator
   {
      public:
         DenseBasisOperator()
         {}

         DenseBasisOperator(int rows, int cols, int dofs)
            :  m_rows(rows),
               m_cols(cols),
               m_dofs(dofs),
               m_data(dofs, mfem::DenseMatrix(rows, cols))
         {
            assert(rows > 0);
            assert(dofs > 0);
            assert(cols > 0);
         }

         DenseBasisOperator(const DenseBasisOperator& other)
            :  BasisOperator(other),
               m_data(other.m_data)
         {}

         DenseBasisOperator(DenseBasisOperator&& other)
            :  BasisOperator(std::move(other)),
               m_data(std::move(other.m_data))
         {}

         DenseBasisOperator& operator=(DenseBasisOperator&& other)
         {
            m_data = std::move(other.m_data);
            return *this;
         }

         void setData(std::vector<mfem::DenseMatrix>&& data)
         {
            m_data = std::move(data);
         }

         mfem::DenseMatrix& operator()(int dof)
         {
            return m_data[dof];
         }

         const mfem::DenseMatrix& operator()(int dof) const
         {
            return m_data[dof];
         }

         int getRows() const override
         {
            return m_rows;
         }

         int getColumns() const override
         {
            return m_cols;
         }

         int getDOFs() const override
         {
            return m_dofs;
         }

         DenseBasisOperator& operator*=(double s) override
         {
            for (int i = 0; i < getDOFs(); i++)
               (*this)(i) *= s;
            return *this;
         }

         double& operator()(int row, int col, int dof)
         {
            return m_data[dof](row, col);
         }

         double operator()(int row, int col, int dof) const override
         {
            return m_data[dof](row, col);
         }

         bool isDense() const override
         {
            return true;
         }

         bool isSparse() const override
         {
            return false;
         }

         std::unique_ptr<BasisOperator>
         MatrixDot(const mfem::DenseMatrix& rhs) const override
         {
            assert(getRows() == rhs.NumRows());
            assert(getColumns() == rhs.NumCols());
            auto result = new DenseBasisOperator(1, 1, getDOFs());
            for (int i = 0; i < getDOFs(); i++)
               (*result)(0, 0, i) = (*this)(i) * rhs;
            return std::unique_ptr<BasisOperator>(result);
         }

         std::unique_ptr<BasisOperator>
         LeftMatrixMult(const mfem::DenseMatrix& lhs) const override
         {
            assert(lhs.NumCols() == getRows());
            auto result = new DenseBasisOperator(lhs.NumRows(), getColumns(), getDOFs());
            for (int i = 0; i < getDOFs(); i++)
               mfem::Mult(lhs, (*this)(i), (*result)(i));
            return std::unique_ptr<BasisOperator>(result);
         }

         std::unique_ptr<BasisOperator>
         RightMatrixMult(const mfem::DenseMatrix& rhs) const override
         {
            assert(getColumns() == rhs.NumRows());
            auto result = new DenseBasisOperator(getRows(), rhs.NumCols(), getDOFs());
            for (int i = 0; i < getDOFs(); i++)
               mfem::Mult((*this)(i), rhs, (*result)(i));
            return std::unique_ptr<BasisOperator>(result);
         }
      private:
         int m_rows;
         int m_cols;
         int m_dofs;
         std::vector<mfem::DenseMatrix> m_data;
   };

   class SparseBasisOperator : public BasisOperator
   {
      public:
         SparseBasisOperator(int rows, int cols, int dofs)
            :  m_rows(rows),
               m_columns(cols),
               m_dofs(dofs)
         {}

         SparseBasisOperator(const SparseBasisOperator& other)
            : BasisOperator(other),
              m_data(other.m_data)
         {}

         SparseBasisOperator(SparseBasisOperator&& other)
            : BasisOperator(std::move(other)),
              m_data(std::move(other.m_data))
         {}

         void setData(std::vector<mfem::SparseMatrix>&& data)
         {
            m_data = std::move(data);
         }

         mfem::SparseMatrix& operator()(int k)
         {
            return m_data[k];
         }

         const mfem::SparseMatrix& operator()(int k) const
         {
            assert(k >= 0);
            assert(k < getDOFs());
            return m_data[k];
         }

         double operator()(int row, int col, int dof) const override
         {
            assert(row >= 0);
            assert(col >= 0);
            assert(dof >= 0);
            assert(row < getRows());
            assert(col < getColumns());
            assert(dof < getDOFs());
            return m_data[dof](row, col);
         }

         SparseBasisOperator& operator*=(double s) override
         {
            for (auto& sm : m_data)
               sm *= s;
            return *this;
         }

         int getRows() const override
         {
            return m_rows;
         }

         int getColumns() const override
         {
            return m_columns;
         }

         int getDOFs() const override
         {
            return m_dofs;
         }

         bool isDense() const override
         {
            return false;
         }

         bool isSparse() const override
         {
            return true;
         }

      private:
         int m_rows;
         int m_columns;
         int m_dofs;
         std::vector<mfem::SparseMatrix> m_data;
   };

   /**
    * @brief Function Scalar Basis Operator (FunctionSBO)
    *
    * Optimized version for functions whose original representation
    * is of the form:
    * @f[
    *    u(x) = \left(
    *       \sum^n_{i=1} w_{1, i} \phi_i(x), \ldots, \sum^n_{i=1} w_{d, i} \phi_i(x) \right)
    * @f]
    */
   class FunctionSBO : public SparseBasisOperator
   {
      public:
         FunctionSBO(mfem::Vector&& shape, int vdim);

         FunctionSBO(const FunctionSBO& other);

         FunctionSBO(FunctionSBO&& other);

         ~FunctionSBO();

      private:
         mfem::Vector m_shape;
         int* m_colIndex;
         std::vector<int*> m_rowIndex;
   };

   /**
    * @brief Jacobian of Scalar Basis Operator (JacobianSBO)
    */
   class JacobianSBO : public SparseBasisOperator
   {
      public:
         /**
          * @brief Constructs the JacobianSBO object.
          * @param[in] dshape Shape gradients in row-major storage order.
          */
         JacobianSBO(mfem::DenseMatrix&& dshape, int sdim, int vdim);

         JacobianSBO(const JacobianSBO& other);

         JacobianSBO(JacobianSBO&& other);

         ~JacobianSBO();

      private:
         mfem::DenseMatrix m_dshape;
         int m_sdim, m_vdim;
         std::vector<int*> m_colIndex;
         int* m_rowIndex;
   };


}

#endif
