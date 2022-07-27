#ifndef RODIN_VARIATIONAL_RANK3OPERATOR_H
#define RODIN_VARIATIONAL_RANK3OPERATOR_H

#include <vector>
#include <mfem.hpp>

namespace Rodin::Variational::Internal
{
   /**
    * @brief Rank-3 tensor for discretized functions on finite bases
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
   class Rank3Operator
   {
      public:
         virtual
         ~Rank3Operator() = default;

         /**
          * @f[
          *    \text{tr} \ A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{p \times p} @f$.
          */
         virtual
         std::unique_ptr<Rank3Operator> Trace() const;

         /**
          * @f[
          *    A(u)^T
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{p \times q} @f$.
          */
         virtual
         std::unique_ptr<Rank3Operator> Transpose() const;

         /**
          * @f$ u @f$ with @f$ n @f$ degrees of freedom
          * @f$ v @f$ with @f$ m @f$ degrees of freedom
          *
          * Stiffness matrix of @f$ m \times n @f$.
          */
         virtual
         std::unique_ptr<Rank3Operator> OperatorSum(const Rank3Operator& rhs) const;

         /**
          * @f[
          *    \Lambda A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         virtual
         std::unique_ptr<Rank3Operator> ScalarMatrixMult(const mfem::DenseMatrix& lhs) const;

         /**
          * @f[
          *    \Lambda A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R} @f$,
          * @f$ \Lambda \in \mathbb{R}^p @f$.
          */
         virtual
         std::unique_ptr<Rank3Operator> ScalarVectorMult(const mfem::Vector& lhs) const;

         /**
          * @f[
          *    A(u) \Lambda
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{s \times p} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         virtual
         std::unique_ptr<Rank3Operator> RightMatrixMult(const mfem::DenseMatrix& rhs) const;

         /**
          * @f[
          *    \Lambda A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{q \times s} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         virtual
         std::unique_ptr<Rank3Operator> LeftMatrixMult(const mfem::DenseMatrix& lhs) const;

         /**
          * @f[
          *    \vec{\lambda} \cdot A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R^d} @f$,
          * @f$ \vec{\lambda} \in \mathbb{R}^{d} @f$.
          */
         virtual
         std::unique_ptr<Rank3Operator> VectorDot(const mfem::Vector& rhs) const;

         /**
          * @f[
          *    \Lambda : A(u)
          * @f]
          * with @f$ A(u) \in \mathbb{R}^{p \times q} @f$,
          * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
          */
         virtual
         std::unique_ptr<Rank3Operator> MatrixDot(const mfem::DenseMatrix& rhs) const;

         /**
          * @f$ u @f$ with @f$ n @f$ degrees of freedom
          * @f$ v @f$ with @f$ m @f$ degrees of freedom
          *
          * Stiffness matrix of @f$ m \times n @f$.
          */
         virtual
         mfem::DenseMatrix OperatorDot(const Rank3Operator& rhs) const;

         /**
          * @brief Adds to element vector.
          */
         virtual
         void AddToVector(mfem::Vector& vec) const;

         virtual
         int GetRows() const = 0;

         virtual
         int GetColumns() const = 0;

         virtual
         int GetDOFs() const = 0;

         virtual
         Rank3Operator& operator=(double s) = 0;

         virtual
         Rank3Operator& operator*=(double s) = 0;

         virtual
         double operator()(int row, int col, int dof) const = 0;
   };

   class DenseRank3Operator : public Rank3Operator
   {
      public:
         DenseRank3Operator()
         {}

         DenseRank3Operator(int rows, int cols, int dofs)
            : m_tensor(rows, cols, dofs)
         {
            assert(rows > 0);
            assert(dofs > 0);
            assert(cols > 0);
         }

         DenseRank3Operator(const DenseRank3Operator& other)
            :  Rank3Operator(other),
               m_tensor(other.m_tensor)
         {}

         DenseRank3Operator(DenseRank3Operator&& other)
            :  Rank3Operator(std::move(other)),
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

         DenseRank3Operator& operator=(double s) override
         {
            for (int i = 0; i < GetDOFs(); i++)
               (*this)(i) = s;
            return *this;
         }

         DenseRank3Operator& operator*=(double s) override
         {
            for (int i = 0; i < GetDOFs(); i++)
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

         std::unique_ptr<Rank3Operator> Trace() const override
         {
            assert(GetRows() == GetColumns());
            auto result = new DenseRank3Operator(1, 1, GetDOFs());
            for (int k = 0; k < GetDOFs(); k++)
               (*result)(0, 0, k) = (*this)(k).Trace();
            return std::unique_ptr<Rank3Operator>(result);
         }

         std::unique_ptr<Rank3Operator>
         MatrixDot(const mfem::DenseMatrix& rhs) const override
         {
            assert(GetRows() == rhs.NumRows());
            assert(GetColumns() == rhs.NumCols());
            auto result = new DenseRank3Operator(1, 1, GetDOFs());
            for (int i = 0; i < GetDOFs(); i++)
               (*result)(0, 0, i) = (*this)(i) * rhs;
            return std::unique_ptr<Rank3Operator>(result);
         }

         std::unique_ptr<Rank3Operator>
         LeftMatrixMult(const mfem::DenseMatrix& lhs) const override
         {
            assert(lhs.NumCols() == GetRows());
            auto result = new DenseRank3Operator(lhs.NumRows(), GetColumns(), GetDOFs());
            for (int i = 0; i < GetDOFs(); i++)
               mfem::Mult(lhs, (*this)(i), (*result)(i));
            return std::unique_ptr<Rank3Operator>(result);
         }

         std::unique_ptr<Rank3Operator>
         RightMatrixMult(const mfem::DenseMatrix& rhs) const override
         {
            assert(GetColumns() == rhs.NumRows());
            auto result = new DenseRank3Operator(GetRows(), rhs.NumCols(), GetDOFs());
            for (int i = 0; i < GetDOFs(); i++)
               mfem::Mult((*this)(i), rhs, (*result)(i));
            return std::unique_ptr<Rank3Operator>(result);
         }
      private:
         mfem::DenseTensor m_tensor;
   };

   /**
    * @brief Optimized version for functions whose original representation
    * is of the form:
    * @f[
    *    u(x) = \left(
    *       \sum^n_{i=1} w_{1, i} \phi_i(x), \ldots, \sum^n_{i=1} w_{d, i} \phi_i(x) \right)
    * @f]
    */
   class ScalarShapeR3O : public Rank3Operator
   {
      public:
         ScalarShapeR3O(mfem::Vector shape, int vdim)
            : m_shape(shape),
              m_vdim(vdim)
         {}

         int GetRows() const override;

         int GetColumns() const override;

         int GetDOFs() const override;

         ScalarShapeR3O& operator*=(double s) override;

         ScalarShapeR3O& operator=(double s) override;

         double operator()(int row, int col, int dof) const override;

         std::unique_ptr<Rank3Operator>
         VectorDot(const mfem::Vector& rhs) const override;
      private:
         mfem::Vector m_shape;
         int m_vdim;
   };

   class JacobianShapeR3O : public Rank3Operator
   {
      public:
         JacobianShapeR3O(mfem::DenseMatrix dshape, int sdim, int vdim)
            : m_dshape(dshape),
              m_sdim(sdim),
              m_vdim(vdim)
         {}

         JacobianShapeR3O(const JacobianShapeR3O& other)
            :  Rank3Operator(other),
               m_dshape(other.m_dshape),
               m_sdim(other.m_sdim),
               m_vdim(other.m_vdim)
         {}

         JacobianShapeR3O(JacobianShapeR3O&& other)
            :  Rank3Operator(std::move(other)),
               m_dshape(std::move(other.m_dshape)),
               m_sdim(other.m_sdim),
               m_vdim(other.m_vdim)
         {}

         int GetRows() const override;

         int GetColumns() const override;

         int GetDOFs() const override;

         JacobianShapeR3O& operator=(double s) override;

         JacobianShapeR3O& operator*=(double s) override;

         double operator()(int row, int col, int dof) const override;

         std::unique_ptr<Rank3Operator> Trace() const override;

      private:
         mfem::DenseMatrix m_dshape;
         int m_sdim, m_vdim;
   };
}

#endif
