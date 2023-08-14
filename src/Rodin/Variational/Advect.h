/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_ADVECT_H
#define RODIN_VARIATIONAL_ADVECT_H

#include <mfem.hpp>

#include "GridFunction.h"
#include "ScalarFunction.h"

namespace Rodin::Variational::Internal
{
  /**
   * @internal
   */
  class DGSolver : public mfem::Solver
  {
    public:
      DGSolver(
          mfem::SparseMatrix& M,
          mfem::SparseMatrix& K,
          const mfem::FiniteElementSpace &fes);

      void SetTimeStep(double dt);

      void SetOperator(const Operator &op);

      virtual void Mult(const mfem::Vector& x, mfem::Vector& y) const;

    private:
      mfem::SparseMatrix& m_M;
      mfem::SparseMatrix& m_K;
      mfem::SparseMatrix m_A;
      mfem::GMRESSolver m_linearSolver;
      mfem::BlockILU m_prec;
      double m_dt;
  };

  /**
   * @internal
   * @brief Time dependent operator for the advection evolution
   *
   * A time-dependent operator for the right-hand side of the ODE. The DG weak
   * form of
   * @f[
   *   \dfrac{du}{dt} = -v \cdot \nabla u
   * @f]
   * is
   * @f[
   *   M \dfrac{du}{dt} = K u + b
   * @f]
   * where @f$ M @f$ and @f$ K @f$ are the
   * mass and advection matrices, and @f$ b @f$ describes the flow on the
   * boundary.
   * This can be written as a general ODE
   * @f[
   * \dfrac{du}{dt} = M^{-1} (K u + b)
   * @f]
   * and this class is used to evaluate the right-hand side.
   */
  class AdvectEvolution : public mfem::TimeDependentOperator
  {
    public:
      AdvectEvolution(mfem::BilinearForm& M,
          mfem::BilinearForm& K, const mfem::Vector& b);

      virtual void Mult(
          const mfem::Vector& x, mfem::Vector& y) const;
      virtual void ImplicitSolve(
          const double dt, const mfem::Vector& x, mfem::Vector& k);

    private:
      mfem::BilinearForm& m_M;
      mfem::BilinearForm& m_K;
      const mfem::Vector& m_b;
      mfem::CGSolver m_MSolver;
      std::unique_ptr<mfem::Solver> m_MPrec;
      std::unique_ptr<DGSolver> m_dgSolver;
  };
}

namespace Rodin::Variational
{
  template <class FES>
  class Advect
  {
    public:
      Advect(GridFunction<FES>& u, const VectorFunctionBase& vel)
        : m_u(u), m_vel(vel.copy())
      {}

      void step(double dt)
      {
        if (!m_initialized)
        {
          m_odeSolver = std::make_unique<mfem::BackwardEulerSolver>();
          m_dgFec    = std::make_unique<mfem::L2_FECollection>(
                      m_u.getFiniteElementSpace().getOrder(),
                      m_u.getFiniteElementSpace().getMesh().getDimension(),
                      mfem::BasisType::GaussLobatto);
          m_dgFes    = std::make_unique<mfem::FiniteElementSpace>(
                      &m_u.getFiniteElementSpace().getMesh().getHandle(),
                      m_dgFec.get());
          m_dgGf    = std::make_unique<mfem::GridFunction>(m_dgFes.get());
          m_bfM     = std::make_unique<mfem::BilinearForm>(m_dgFes.get());
          m_bfK     = std::make_unique<mfem::BilinearForm>(m_dgFes.get());
          m_lf      = std::make_unique<mfem::LinearForm>(m_dgFes.get());
          m_solution  = std::make_unique<mfem::GridFunction>(m_dgFes.get());
          m_uCoeff   = std::make_unique<mfem::GridFunctionCoefficient>(&m_u.getHandle());
          m_velCoeff  = m_vel->build();
          std::cout << m_velCoeff->GetVDim() << std::endl;
          m_solutionCoeff  = std::make_unique<mfem::GridFunctionCoefficient>(m_solution.get());
          m_transfer      = std::make_unique<mfem::TransferOperator>(
                          m_u.getFiniteElementSpace().getFES(), *m_dgFes);

          // Build M
          m_bfM->AddDomainIntegrator(new mfem::MassIntegrator);

          // Build K
          constexpr double alpha = -1.0;
          m_bfK->AddDomainIntegrator(
            new mfem::ConvectionIntegrator(*m_velCoeff, alpha));
          m_bfK->AddInteriorFaceIntegrator(
            new mfem::NonconservativeDGTraceIntegrator(*m_velCoeff, alpha));
          m_bfK->AddBdrFaceIntegrator(
            new mfem::NonconservativeDGTraceIntegrator(*m_velCoeff, alpha));

          // Build linear form
          mfem::ConstantCoefficient inflow(0.0);
          m_lf->AddBdrFaceIntegrator(
              new mfem::BoundaryFlowIntegrator(inflow, *m_velCoeff, alpha));

          // Assemble all the forms
          m_bfM->Assemble();
          bool skipZeros = false;
          m_bfK->Assemble(skipZeros);
          m_lf->Assemble();
          m_bfM->Finalize();
          m_bfK->Finalize(skipZeros);

          // Create the advect evolution
          m_adv = std::make_unique<Internal::AdvectEvolution>(*m_bfM, *m_bfK, *m_lf);
          m_adv->SetTime(0.0);
          m_odeSolver->Init(*m_adv);

          m_initialized = true;
        }

        // Transfer u into the L2 space
        m_solution->ProjectCoefficient(*m_uCoeff);

        // Make the step
        m_odeSolver->Step(*m_solution, m_t, dt);

        // Transfer back
        m_u.getHandle().ProjectCoefficient(*m_solutionCoeff);

        m_t += dt;
      }

      bool isInitialized() const
      {
        return m_initialized;
      }

      const GridFunction<FES>& getSolution() const
      {
        return m_u;
      }

      const GridFunction<FES>& getVelocity() const
      {
        return m_vel;
      }

    private:
      double m_t;
      GridFunction<FES>& m_u;
      std::unique_ptr<VectorFunctionBase> m_vel;

      bool m_initialized;
      std::unique_ptr<mfem::ODESolver> m_odeSolver;
      std::unique_ptr<mfem::L2_FECollection> m_dgFec;
      std::unique_ptr<mfem::FiniteElementSpace> m_dgFes;
      std::unique_ptr<mfem::GridFunction> m_dgGf;
      std::unique_ptr<mfem::BilinearForm> m_bfM;
      std::unique_ptr<mfem::BilinearForm> m_bfK;
      std::unique_ptr<mfem::LinearForm>  m_lf;
      std::unique_ptr<mfem::GridFunction> m_solution;
      std::unique_ptr<mfem::GridFunctionCoefficient> m_uCoeff;
      std::unique_ptr<mfem::VectorCoefficient> m_velCoeff;
      std::unique_ptr<mfem::GridFunctionCoefficient> m_solutionCoeff;
      std::unique_ptr<mfem::TransferOperator> m_transfer;
      std::unique_ptr<Internal::AdvectEvolution> m_adv;
  };
  template <class FES>
  Advect(GridFunction<FES>&, const VectorFunctionBase&) -> Advect<FES>;
}

#endif
