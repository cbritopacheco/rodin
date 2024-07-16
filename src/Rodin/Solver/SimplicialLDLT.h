/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SIMPLICIALLDLT_H
#define RODIN_SOLVER_SIMPLICIALLDLT_H

#include <Eigen/SparseCholesky>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  /**
   * @defgroup SimplicialLDLTSpecializations SimplicialLDLT Template Specializations
   * @brief Template specializations of the SimplicialLDLT class.
   * @see SimplicialLDLT
   */

  /**
   * @ingroup SimplicialLDLTSpecializations
   * @brief A direct sparse LDLT Cholesky factorizations without square root
   * for use with Math::SparseMatrix<Real> and Math::Vector<Real>.
   */
  template <class Scalar>
  class SimplicialLDLT<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>> final
    : public SolverBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using ProblemType = Variational::ProblemBase<OperatorType, VectorType, ScalarType>;

      using Parent = SolverBase<OperatorType, VectorType, ScalarType>;

      using Parent::solve;

      SimplicialLDLT(ProblemType& pb)
        : Parent(pb)
      {}

      SimplicialLDLT(const SimplicialLDLT& other)
        : Parent(other)
      {}

      SimplicialLDLT(SimplicialLDLT&& other)
        : Parent(std::move(other))
      {}

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        x = m_solver.compute(A).solve(b);
      }

      inline
      SimplicialLDLT* copy() const noexcept override
      {
        return new SimplicialLDLT(*this);
      }

    private:
      Eigen::SimplicialLDLT<OperatorType> m_solver;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for SimplicialLDLT
   */
  template <class Scalar>
  SimplicialLDLT(Variational::ProblemBase<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>, Scalar>&)
    -> SimplicialLDLT<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>;

}

#endif


