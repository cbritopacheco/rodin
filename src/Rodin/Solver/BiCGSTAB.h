/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_BiCGSTAB_H
#define RODIN_SOLVER_BiCGSTAB_H

#include <optional>
#include <functional>
#include <Eigen/SparseCholesky>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup BiCGSTABSpecializations BiCGSTAB Template Specializations
   * @brief Template specializations of the BiCGSTAB class.
   * @see BiCGSTAB
   */

  /**
   * @ingroup BiCGSTABSpecializations
   * @brief Conjugate gradient solver for self-adjoint problems, for use with
   * Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class BiCGSTAB<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      /**
       * @brief Constructs the BiCGSTAB object with default parameters.
       */
      BiCGSTAB(ProblemType& pb)
        : Parent(pb)
      {}

      BiCGSTAB(const BiCGSTAB& other)
        : Parent(other)
      {}

      BiCGSTAB(BiCGSTAB&& other)
        : Parent(std::move(other))
      {}

      ~BiCGSTAB() = default;

      BiCGSTAB& setTolerance(double tol)
      {
        m_solver.setTolerance(tol);
        return *this;
      }

      BiCGSTAB& setMaxIterations(size_t maxIt)
      {
        m_solver.setMaxIterations(maxIt);
        return *this;
      }

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        m_solver.compute(A);
        x = m_solver.solve(b);
      }

      bool success() const
      {
        return m_solver.info() == Eigen::Success;
      }

      BiCGSTAB* copy() const noexcept override
      {
        return new BiCGSTAB(*this);
      }

    private:
      Eigen::BiCGSTAB<OperatorType> m_solver;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for BiCGSTAB
   */
  template <class Scalar>
  BiCGSTAB(Variational::ProblemBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> BiCGSTAB<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;
}

#endif


