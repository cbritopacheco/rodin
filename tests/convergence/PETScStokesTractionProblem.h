/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_STOKESTRACTIONPROBLEM_H
#define RODIN_TESTS_CONVERGENCE_PETSC_STOKESTRACTIONPROBLEM_H

#include <functional>
#include <utility>
#include <gtest/gtest.h>

#include "Stokes.h"
#include "Rodin/PETSc.h"
#ifdef RODIN_USE_MPI
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif

namespace Rodin::Tests::Convergence
{
  struct StokesTractionErrors
  {
      StokesErrors fields;
      Real pressureMean;
  };

  /** @brief Physical-stress Stokes verification with mixed boundary data.
   * @par Mathematical formulation
   * @f$-\nabla\cdot(2\varepsilon(u)-pI)=f@f$ and
   * @f$\nabla\cdot u=0@f$, with velocity on the Dirichlet partition and
   * @f$(2\varepsilon(u)-pI)n@f$ on the traction partition.
   * Prescribed traction fixes the pressure level: no pressure-mean multiplier
   * or nullspace constraint is introduced. A constant shift of the prescribed
   * traction by @f$c n@f$ changes the solution pressure to @f$p-c@f$ without
   * changing velocity or the volume source.
   * @par Architecture
   * Each solve owns fresh Taylor--Hood spaces and a two-field PETSc system.
   * MUMPS factorization, convergence status, an independent coefficient
   * residual, physical error norms, and the pressure integral are checked or
   * measured separately. Pointwise stress evaluation is noncollective.
   */
  template <class MeshType>
  class PETScStokesTractionProblem
  {
    public:
      static constexpr Geometry::Attribute DirichletAttribute = 401;
      static constexpr Geometry::Attribute TractionAttribute = 402;

      PETScStokesTractionProblem(
        const MeshType& mesh, const StokesData& data, size_t assemblyOrder = 16)
        : m_mesh(mesh),
          m_data(data),
          m_order(assemblyOrder)
      {}

      /** @brief Solve and measure the global mixed-boundary problem.
       * @note All ranks of an MPI mesh communicator must participate in the
       * assembly, solve, residual norms, field norms, and pressure integral.
       * The traction pressure shift changes boundary data only; it is a
       * negative-control parameter, not an alternative pressure constraint.
       */
      template <size_t K>
      StokesTractionErrors solve(
        size_t normOrder = 18, Real tractionPressureShift = 0) const
      {
        static_assert(K >= 2);
        using namespace Variational;
        const auto& mesh = m_mesh.get();
        H1<K, Math::SpatialVector<Real>, MeshType> velocitySpace(
          std::integral_constant<size_t, K>{}, mesh, mesh.getDimension());
        H1<K - 1, Real, MeshType> pressureSpace(
          std::integral_constant<size_t, K - 1>{}, mesh);
        PETSc::Variational::TrialFunction u(velocitySpace);
        PETSc::Variational::TrialFunction p(pressureSpace);
        PETSc::Variational::TestFunction v(velocitySpace);
        PETSc::Variational::TestFunction q(pressureSpace);
        const BoundaryNormal normal(mesh);
        const VectorFunction traction(mesh.getDimension(),
          [data = m_data, normal, tractionPressureShift](const Geometry::Point& point) {
            const auto stress = data.getStress(point.getPhysicalCoordinates());
            const auto outward = normal(point);
            Math::SpatialVector<Real> value(static_cast<std::uint8_t>(outward.size()));
            value.setZero();
            for (size_t i = 0; i < value.size(); ++i)
            {
              for (size_t j = 0; j < value.size(); ++j)
                value(i) += stress(i, j) * outward(j);
              value(i) += tractionPressureShift * outward(i);
            }
            return value;
          });
        auto diffusion = Integral(
          Jacobian(u) + Jacobian(u).T(), Real(0.5) * (Jacobian(v) + Jacobian(v).T()));
        auto pressureVelocity = Integral(p, Div(v));
        auto incompressibility = Integral(Div(u), q);
        auto body = Integral(m_data.getForcing(), v);
        auto boundaryLoad = BoundaryIntegral(traction, v);
        diffusion.setOrder(m_order);
        pressureVelocity.setOrder(m_order);
        incompressibility.setOrder(m_order);
        body.setOrder(m_order);
        boundaryLoad.setOrder(m_order);
        Problem problem(u, p, v, q);
        problem = diffusion - pressureVelocity + incompressibility - body -
          boundaryLoad.over(TractionAttribute) +
          DirichletBC(u, m_data.getVelocity()).on(DirichletAttribute);
        for (const auto& [name, value] : {std::pair{"-mat_mumps_icntl_14", "100"},
               std::pair{"-mat_mumps_icntl_20", "0"}})
        {
          PetscBool set = PETSC_FALSE;
          EXPECT_EQ(PetscOptionsHasName(nullptr, nullptr, name, &set), PETSC_SUCCESS);
          if (!set)
          {
            EXPECT_EQ(PetscOptionsSetValue(nullptr, name, value), PETSC_SUCCESS);
          }
        }
        PETSc::Solver::KSP solver(problem);
        solver.setType(KSPPREONLY);
        PC pc = nullptr;
        EXPECT_EQ(KSPGetPC(solver.getHandle(), &pc), PETSC_SUCCESS);
        EXPECT_EQ(PCSetType(pc, PCLU), PETSC_SUCCESS);
        EXPECT_EQ(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS), PETSC_SUCCESS);
        solver.solve();
        PCFailedReason failure = PC_NOERROR;
        EXPECT_EQ(PCGetFailedReason(pc, &failure), PETSC_SUCCESS);
        EXPECT_EQ(failure, PC_NOERROR);
        Mat factor = nullptr;
        EXPECT_EQ(PCFactorGetMatrix(pc, &factor), PETSC_SUCCESS);
        PetscInt info = 0, detail = 0;
        EXPECT_EQ(MatMumpsGetInfog(factor, 1, &info), PETSC_SUCCESS);
        EXPECT_EQ(MatMumpsGetInfog(factor, 2, &detail), PETSC_SUCCESS);
        EXPECT_EQ(info, 0) << "MUMPS INFOG(2)=" << detail;
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
        EXPECT_EQ(KSPGetConvergedReason(solver.getHandle(), &reason), PETSC_SUCCESS);
        EXPECT_GT(reason, 0);
        const auto& system = problem.getLinearSystem();
        Vec residual = nullptr;
        EXPECT_EQ(VecDuplicate(system.getVector(), &residual), PETSC_SUCCESS);
        EXPECT_EQ(
          MatMult(system.getOperator(), system.getSolution(), residual), PETSC_SUCCESS);
        EXPECT_EQ(VecAXPY(residual, -1, system.getVector()), PETSC_SUCCESS);
        PetscReal norm = 0, rhsNorm = 0;
        EXPECT_EQ(VecNorm(residual, NORM_2, &norm), PETSC_SUCCESS);
        EXPECT_EQ(VecNorm(system.getVector(), NORM_2, &rhsNorm), PETSC_SUCCESS);
        EXPECT_TRUE(std::isfinite(norm));
        EXPECT_TRUE(std::isfinite(rhsNorm));
        EXPECT_LT(norm / std::max(Real(1), rhsNorm), IndependentResidualTolerance);
        EXPECT_EQ(VecDestroy(&residual), PETSC_SUCCESS);
        const auto& velocity = u.getSolution();
        const auto& pressure = p.getSolution();
        auto pressureIntegral = Integral(pressure);
        pressureIntegral.setOrder(normOrder);
        const Real pressureMean = pressureIntegral.compute();
        EXPECT_TRUE(std::isfinite(pressureMean));
        return {{ErrorNorm::computeVector(mesh, velocity, m_data.getVelocity(),
                   m_data.getVelocityJacobian(), normOrder),
                  ErrorNorm::compute(mesh, pressure, m_data.getPressure(),
                    m_data.getPressureGradient(), normOrder),
                  ErrorNorm::computeDivergenceL2(mesh, velocity, normOrder)},
          pressureMean};
      }

    private:
      // Dimensionless coefficient residual; below all field acceptance budgets.
      static constexpr Real IndependentResidualTolerance = 1e-11;
      std::reference_wrapper<const MeshType> m_mesh;
      StokesData m_data;
      size_t m_order;
  };
}

#endif
