/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_UMFPACK_H
#define RODIN_SOLVER_UMFPACK_H

#ifdef RODIN_USE_UMFPACK

#include <optional>
#include <functional>

#include <Eigen/UmfPackSupport>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"

namespace Rodin::Solver
{
  /**
   * @defgroup UMFPackSpecializations UMFPack Template Specializations
   * @brief Template specializations of the UMFPack class.
   * @see UMFPack
   */

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for UMFPack
   */
  template <class LinearSystem>
  UMFPack(Variational::ProblemBase<LinearSystem>&) -> UMFPack<LinearSystem>;

  /**
   * @ingroup UMFPackSpecializations
   * @brief UMFPack for use with Math::SparseMatrix and Math::Vector.
   */
  template <class Scalar>
  class UMFPack<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
    : public SolverBase<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>
  {
    public:
      using ScalarType = Scalar;

      using VectorType = Math::Vector<ScalarType>;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using LinearSystemType = Math::LinearSystem<OperatorType, VectorType>;

      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;

      using Parent = SolverBase<LinearSystemType>;

      using Parent::solve;

      UMFPack(ProblemBaseType& pb)
        : Parent(pb)
      {}

      UMFPack(const UMFPack& other)
        : Parent(other)
      {}

      UMFPack(UMFPack&& other)
        : Parent(std::move(other))
      {}

      ~UMFPack() = default;

      void solve(LinearSystemType& axb) override
      {
        m_solver.compute(axb.getOperator());
        axb.getSolution() = m_solver.solve(axb.getVector());
      }

      void printControl()
      {
        m_solver.printUmfpackControl();
      }

      void printInfo()
      {
        m_solver.printUmfpackInfo();
      }

      void printStatus()
      {
        m_solver.printUmfpackStatus();
      }

      UMFPack* copy() const noexcept override
      {
        return new UMFPack(*this);
      }

    private:
      Eigen::UmfPackLU<OperatorType> m_solver;
  };
}

#endif // #ifdef RODIN_USE_UMFPACK
#endif


