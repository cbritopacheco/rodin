/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Contains the implementation for the Advect class
 *
 * The internal code for the Advect functionality has been taken and adapted
 * from the mfem library, example 9. This code is distributed under the BSD
 * 3-Clause "New" or "Revised" License.
 *
 * @see https://github.com/mfem/mfem/blob/master/LICENSE
 * @see https://github.com/mfem/mfem/blob/master/examples/ex9.cpp
 */
#include "Advect.h"

namespace Rodin::Variational::Internal
{
  DGSolver::DGSolver(
      mfem::SparseMatrix& M,
      mfem::SparseMatrix& K,
      const mfem::FiniteElementSpace &fes)
    : m_M(M), m_K(K),
      m_prec(fes.GetFE(0)->GetDof(),
         mfem::BlockILU::Reordering::MINIMUM_DISCARDED_FILL),
      m_dt(-1.0)
  {
    m_linearSolver.iterative_mode = false;
    m_linearSolver.SetRelTol(1e-9);
    m_linearSolver.SetAbsTol(0.0);
    m_linearSolver.SetMaxIter(100);
    m_linearSolver.SetPrintLevel(0);
    m_linearSolver.SetPreconditioner(m_prec);
  }

  void DGSolver::SetTimeStep(double dt)
  {
    if (dt != m_dt)
    {
      m_dt = dt;

      // Form operator A = M - dt * K
      m_A = m_K;
      m_A *= -m_dt;
      m_A += m_M;

      // This will also call SetOperator on the preconditioner
      m_linearSolver.SetOperator(m_A);
    }
  }

  void DGSolver::SetOperator(const Operator &op)
  {
    m_linearSolver.SetOperator(op);
  }

  void DGSolver::Mult(const mfem::Vector& x, mfem::Vector& y) const
  {
    m_linearSolver.Mult(x, y);
  }

  AdvectEvolution::AdvectEvolution(
      mfem::BilinearForm& M,
      mfem::BilinearForm& K,
      const mfem::Vector& b)
    : TimeDependentOperator(M.Height()), m_M(M), m_K(K), m_b(b)
  {
    mfem::Array<int> ess_tdof_list;
    if (M.GetAssemblyLevel() == mfem::AssemblyLevel::LEGACY)
    {
      m_MPrec = std::make_unique<mfem::DSmoother>(M.SpMat());
      m_MSolver.SetOperator(M.SpMat());
      m_dgSolver = std::make_unique<DGSolver>(M.SpMat(), K.SpMat(), *M.FESpace());
    }
    else
    {
      m_MPrec = std::make_unique<mfem::OperatorJacobiSmoother>(M, ess_tdof_list);
      m_MSolver.SetOperator(M);
      m_dgSolver = nullptr;
    }
    m_MSolver.SetPreconditioner(*m_MPrec);
    m_MSolver.iterative_mode = false;
    m_MSolver.SetRelTol(1e-9);
    m_MSolver.SetAbsTol(0.0);
    m_MSolver.SetMaxIter(100);
    m_MSolver.SetPrintLevel(0);
  }

  void AdvectEvolution::Mult(const mfem::Vector& x, mfem::Vector& y) const
  {
    // y = M^{-1} (K x + b)
    mfem::Vector z(m_M.Height());
    m_K.Mult(x, z);
    z += m_b;
    m_MSolver.Mult(z, y);
  }

  void AdvectEvolution::ImplicitSolve(const double dt,
      const mfem::Vector &x, mfem::Vector &k)
  {
    MFEM_VERIFY(m_dgSolver != NULL,
            "Implicit time integration is not supported with partial assembly");

    mfem::Vector z(m_M.Height());
    m_K.Mult(x, z);
    z += m_b;
    m_dgSolver->SetTimeStep(dt);
    m_dgSolver->Mult(z, k);
  }
}
