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
         /**
          * @f[
          *    \text{tr} \ A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{p \times p} @f$.
          */
         virtual
         std::unique_ptr<BasisOperator> Trace() const;

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

         /**
          * @f$ u @f$ with @f$ n @f$ degrees of freedom
          * @f$ v @f$ with @f$ m @f$ degrees of freedom
          *
          * Stiffness matrix of @f$ m \times n @f$.
          */
         virtual
         mfem::DenseMatrix OperatorDot(const BasisOperator& rhs) const;

         constexpr bool isDense() const
         {
            return !isSparse();
         }

         virtual ~BasisOperator() = default;

         /**
          * @brief Adds to element vector.
          */
         virtual void addToVector(mfem::Vector& vec) const;

         virtual bool isSparse() const = 0;

         virtual int getRows() const = 0;

         virtual int getColumns() const = 0;

         virtual int getDOFs() const = 0;

         virtual BasisOperator& operator=(double s) = 0;

         virtual BasisOperator& operator*=(double s) = 0;

         virtual double operator()(int row, int col, int dof) const = 0;
   };

   class DenseBasisOperator : public BasisOperator
   {
      public:
         DenseBasisOperator()
         {}

         DenseBasisOperator(int rows, int cols, int dofs)
            : m_tensor(rows, cols, dofs)
         {
            assert(rows > 0);
            assert(dofs > 0);
            assert(cols > 0);
         }

         DenseBasisOperator(const DenseBasisOperator& other)
            :  BasisOperator(other),
               m_tensor(other.m_tensor)
         {}

         DenseBasisOperator(DenseBasisOperator&& other)
            :  BasisOperator(std::move(other)),
               m_tensor(std::move(other.m_tensor))
         {}

         mfem::DenseMatrix& operator()(int dof)
         {
            return m_tensor(dof);
         }

         const mfem::DenseMatrix& operator()(int dof) const
         {
            return m_tensor(dof);
         }

         int getRows() const override
         {
            return m_tensor.SizeI();
         }

         int getColumns() const override
         {
            return m_tensor.SizeJ();
         }

         int getDOFs() const override
         {
            return m_tensor.SizeK();
         }

         DenseBasisOperator& operator=(double s) override
         {
            for (int i = 0; i < getDOFs(); i++)
               (*this)(i) = s;
            return *this;
         }

         DenseBasisOperator& operator*=(double s) override
         {
            for (int i = 0; i < getDOFs(); i++)
               (*this)(i) *= s;
            return *this;
         }

         double& operator()(int row, int col, int dof)
         {
            return m_tensor(row, col, dof);
         }

         double operator()(int row, int col, int dof) const override
         {
            return m_tensor(row, col, dof);
         }

         bool isSparse() const override
         {
            return false;
         }

         std::unique_ptr<BasisOperator> Trace() const override
         {
            assert(getRows() == getColumns());
            auto result = new DenseBasisOperator(1, 1, getDOFs());
            for (int k = 0; k < getDOFs(); k++)
               (*result)(0, 0, k) = (*this)(k).Trace();
            return std::unique_ptr<BasisOperator>(result);
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
         mfem::DenseTensor m_tensor;
   };

   // class SparseBasisOperator : public BasisOperator
   // {
   //    public:
   //       class SparsityPattern
   //       {
   //          public:
   //             SparsityPattern(
   //                   std::set<std::tuple<int>> 
   //                   int rows, int cols, int dofs);

   //          private:
   //             mfem::Memory<int> m_i;
   //             mfem::Memory<int> m_j;
   //             mfem::Memory<int> m_k;
   //       };

   //       SparseBasisOperator(std::unique_ptr<double[]> data, const SparsityPattern& pattern);

   //       SparsityPattern& getSparsityPattern()
   //       {
   //          return m_sparsityPattern;
   //       }

   //       const SparsityPattern& getSparsityPattern() const
   //       {
   //          return m_sparsityPattern;
   //       }

   //       bool isSparse() const override
   //       {
   //          return true;
   //       }

   //    private:
   //       SparsityPattern m_sparsityPattern;
   //       mfem::Memory<double> m_data;
   // };

   /**
    * @brief Scalar Shape Function Basis Operator (SSFBO)
    *
    * Optimized version for functions whose original representation
    * is of the form:
    * @f[
    *    u(x) = \left(
    *       \sum^n_{i=1} w_{1, i} \phi_i(x), \ldots, \sum^n_{i=1} w_{d, i} \phi_i(x) \right)
    * @f]
    */
   class SSFBO : public BasisOperator
   {
      public:
         SSFBO(const mfem::Vector& shape, int vdim)
            : m_shape(shape),
              m_vdim(vdim)
         {}

         SSFBO(mfem::Vector&& shape, int vdim)
            : m_shape(std::move(shape)),
              m_vdim(vdim)
         {}

         SSFBO(const SSFBO& other)
            :  BasisOperator(other),
               m_shape(other.m_shape),
               m_vdim(other.m_vdim)
         {}

         SSFBO(SSFBO&& other)
            :  BasisOperator(std::move(other)),
               m_shape(std::move(other.m_shape)),
               m_vdim(other.m_vdim)
         {
            other.m_vdim = 0;
         }

         int getRows() const override;

         int getColumns() const override;

         int getDOFs() const override;

         SSFBO& operator*=(double s) override;

         SSFBO& operator=(double s) override;

         double operator()(int row, int col, int dof) const override;

         bool isSparse() const override
         {
            return true;
         }

      private:
         mfem::Vector m_shape;
         int m_vdim;
   };

   /**
    * @brief Jacobian of Scalar Shape Function Basis Operator (JSSFBO)
    */
   class JSSFBO : public BasisOperator
   {
      public:
         JSSFBO(const mfem::DenseMatrix& dshape, int sdim, int vdim)
            : m_dshape(dshape),
              m_sdim(sdim),
              m_vdim(vdim)
         {}

         JSSFBO(mfem::DenseMatrix&& dshape, int sdim, int vdim)
            : m_sdim(sdim),
              m_vdim(vdim)
         {
            m_dshape.SetSize(dshape.Height(), dshape.Width());
            m_dshape.GetMemory() = std::move(dshape.GetMemory());
         }

         JSSFBO(const JSSFBO& other)
            :  BasisOperator(other),
               m_dshape(other.m_dshape),
               m_sdim(other.m_sdim),
               m_vdim(other.m_vdim)
         {}

         JSSFBO(JSSFBO&& other)
            :  BasisOperator(std::move(other)),
               m_sdim(other.m_sdim),
               m_vdim(other.m_vdim)
         {
            m_dshape.SetSize(other.m_dshape.Height(), other.m_dshape.Width());
            m_dshape.GetMemory() = std::move(other.m_dshape.GetMemory());
            other.m_sdim = 0;
            other.m_vdim = 0;
         }

         int getRows() const override;

         int getColumns() const override;

         int getDOFs() const override;

         JSSFBO& operator=(double s) override;

         JSSFBO& operator*=(double s) override;

         double operator()(int row, int col, int dof) const override;

         bool isSparse() const override
         {
            return true;
         }

         std::unique_ptr<BasisOperator> Trace() const override;

      private:
         mfem::DenseMatrix m_dshape;
         int m_sdim, m_vdim;
   };


}

#endif
