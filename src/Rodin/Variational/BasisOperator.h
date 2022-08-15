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

         virtual ~BasisOperator() = default;

         /**
          * @brief Adds to element vector.
          */
         virtual void addToVector(mfem::Vector& vec) const = 0;

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
               m_rows(other.m_rows),
               m_cols(other.m_cols),
               m_dofs(other.m_dofs),
               m_data(other.m_data)
         {}

         DenseBasisOperator(DenseBasisOperator&& other)
            :  BasisOperator(std::move(other)),
               m_rows(other.m_rows),
               m_cols(other.m_cols),
               m_dofs(other.m_dofs),
               m_data(std::move(other.m_data))
         {
            other.m_rows = 0;
            other.m_cols = 0;
            other.m_dofs = 0;
         }

         DenseBasisOperator& operator=(DenseBasisOperator&& other)
         {
            m_data = std::move(other.m_data);
            m_rows = other.m_rows;
            m_cols = other.m_cols;
            m_dofs = other.m_dofs;

            other.m_rows = 0;
            other.m_cols = 0;
            other.m_dofs = 0;
            return *this;
         }

         void transpose()
         {
            std::swap(m_rows, m_cols);
            for (auto& v : m_data)
               v.Transpose();
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

         DenseBasisOperator& operator+=(const DenseBasisOperator& rhs)
         {
            assert(getRows() == rhs.getRows());
            assert(getColumns() == rhs.getColumns());
            assert(getDOFs() == rhs.getDOFs());
            for (int i = 0; i < getDOFs(); i++)
               operator()(i) += rhs(i);
            return *this;
         }

         DenseBasisOperator& operator*=(double s) override
         {
            for (int i = 0; i < getDOFs(); i++)
               operator()(i) *= s;
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

         void addToVector(mfem::Vector& vec) const override;

      private:
         int m_rows;
         int m_cols;
         int m_dofs;
         std::vector<mfem::DenseMatrix> m_data;
   };
}

#endif
