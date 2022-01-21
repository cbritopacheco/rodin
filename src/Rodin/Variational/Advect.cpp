/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Advect.h"
#include "GridFunction.h"
#include "H1.h"
#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
   namespace Internal
   {
      DGSolver::DGSolver(
                  mfem::SparseMatrix& M,
                  mfem::SparseMatrix& K,
                  const mfem::FiniteElementSpace& fes)
         :  m_M(M),
            m_K(K),
            m_prec(
                  fes.GetFE(0)->GetDof(),
                  mfem::BlockILU::Reordering::MINIMUM_DISCARDED_FILL),
            m_dt(-1.0)
      {
         m_linearSolver.iterative_mode = false;
         m_linearSolver.SetRelTol(1e-9);
         m_linearSolver.SetAbsTol(0.0);
         m_linearSolver.SetMaxIter(100);
         m_linearSolver.SetPrintLevel(mfem::IterativeSolver::PrintLevel().Errors());
         m_linearSolver.SetPreconditioner(m_prec);
      }

      void DGSolver::SetTimeStep(double dt)
      {
         if (m_dt != dt)
         {
            mfem::SparseMatrix A;
            m_dt = dt;
            A = m_K;
            A *= -dt;
            A += m_M;

            // This will also call SetOperator on the preconditioner
            m_linearSolver.SetOperator(A);
         }
      }

      void DGSolver::SetOperator(const mfem::Operator& op)
      {
         m_linearSolver.SetOperator(op);
      }

      void DGSolver::Mult(const mfem::Vector& x, mfem::Vector& y) const
      {
         m_linearSolver.Mult(x, y);
      }

      AdvectionEvolution::AdvectionEvolution(
            mfem::BilinearForm& M, mfem::BilinearForm& K)
         :  mfem::TimeDependentOperator(M.Height()),
            m_M(M),
            m_K(K),
            m_z(M.Height())
      {
         m_solver.SetPreconditioner(*m_prec);
         m_solver.iterative_mode = false;
         m_solver.SetRelTol(1e-9);
         m_solver.SetAbsTol(0.0);
         m_solver.SetMaxIter(100);
         m_solver.SetPrintLevel(mfem::IterativeSolver::PrintLevel().Errors());
      }

      void AdvectionEvolution::Mult(
            const mfem::Vector& x, mfem::Vector& y) const
      {
         // y = M^{-1} (K x + b)
         m_K.Mult(x, m_z);
         m_solver.Mult(m_z, y);
      }

      void AdvectionEvolution::ImplicitSolve(
            double dt, const mfem::Vector& x, mfem::Vector& k)
      {
         MFEM_VERIFY(m_dgSolver != nullptr,
               "Implicit time integration is not supported with partial assembly");
         m_K.Mult(x, m_z);
         m_dgSolver->SetTimeStep(dt);
         m_dgSolver->Mult(m_z, k);
      }
   }

   Advect::Advect(GridFunction<H1>& u, const VectorCoefficientBase& velocity)
      :  m_t(0.0),
         m_u(u),
         m_velocity(velocity.copy()),
         m_odeSolver(new mfem::SDIRK33Solver),
         m_fec(3,
               u.getFiniteElementSpace().getMesh().getDimension(),
               mfem::BasisType::GaussLobatto),
         m_fes(&u.getFiniteElementSpace().getMesh().getHandle(), &m_fec),
         m_M(&m_fes),
         m_K(&m_fes),
         m_adv(m_M, m_K)
   {
      assert(u.getFiniteElementSpace().getRangeDimension() == 1);
      m_velocity->buildMFEMVectorCoefficient();
      m_M.AddDomainIntegrator(
            new mfem::ConvectionIntegrator(
               m_velocity->getMFEMVectorCoefficient(), -1.0));
      m_K.AddInteriorFaceIntegrator(
            new mfem::NonconservativeDGTraceIntegrator(
               m_velocity->getMFEMVectorCoefficient(), -1.0));
      m_K.AddBdrFaceIntegrator(
            new mfem::NonconservativeDGTraceIntegrator(
               m_velocity->getMFEMVectorCoefficient(), -1.0));

      m_M.Assemble();
      m_K.Assemble(0);
      m_M.Finalize();
      m_K.Finalize(0);

      m_adv.SetTime(m_t);
      m_odeSolver->Init(m_adv);
   }

   void Advect::step(double dt)
   {
      m_odeSolver->Step(m_u.getHandle(), m_t, dt);
   }

   double Advect::getTime() const
   {
      return m_t;
   }
}
