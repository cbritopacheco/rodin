/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_NONLINEAR_POISSON_H
#define RODIN_TESTS_CONVERGENCE_PETSC_NONLINEAR_POISSON_H

#include <utility>
#include <functional>
#include <variant>
#ifdef RODIN_USE_MPI
#include <boost/mpi/collectives.hpp>
#endif

#include "NonlinearPoisson.h"
#include "Rodin/PETSc.h"

namespace Rodin::Tests::Convergence
{
  /** @brief PETSc/SNES verification workload for semilinear Poisson.
   * @par Architecture
   * One mesh and degree define fixed spaces and PETSc layouts. The continuous
   * manufactured data is shared with the native Newton suites. A separate
   * state field is synchronized from each SNES iterate before assembly;
   * residual and Jacobian callbacks operate on the same variational problem.
   * Each new workload owns new fields and solvers, without resizing a system.
   * The derivative oracle uses SNES callbacks on distinct perturbed vectors.
   * With opt-in lifting, SNES stores the homogeneous correction @f$w_h@f$.
   * Callback updates reconstruct @f$u_h=I_hu+w_h@f$ through grid-function
   * accumulation, including MPI ghost synchronization. Default flat-mesh
   * behavior retains zero lifting.
   */
  template <size_t K, class MeshType>
  class PETScNonlinearPoissonProblem
  {
    public:
      enum class Boundary
      {
        Dirichlet,
        MixedNeumann,
        Robin,
        PureNeumann
      };
      static constexpr Geometry::Attribute DirichletAttribute = 601;
      static constexpr Geometry::Attribute NaturalAttribute = 602;

      using SpaceType = Variational::H1<K, Real, MeshType>;
      using StateType = PETSc::Variational::GridFunction<SpaceType>;
      using TrialType = PETSc::Variational::TrialFunction<StateType, SpaceType>;
      using TestType = PETSc::Variational::TestFunction<SpaceType>;
      using ProblemType =
        Variational::Problem<PETSc::Math::LinearSystem, TrialType, TestType>;

      /** @brief Construct and assemble a fixed-layout nonlinear workload.
       * @note MPI space construction, essential-DOF assembly, global free-DOF
       * counting, and distributed problem assembly require communicator
       * participation. Pointwise manufactured data and normal contractions
       * remain noncollective. Pure Neumann uses the cached global space size;
       * the positive reaction removes the constant nullspace.
       */
      explicit PETScNonlinearPoissonProblem(const MeshType& mesh,
        size_t quadratureOrder = 12, Real amplitude = 1, bool omitCubic = false,
        bool wrongTangent = false, bool liftBoundary = false,
        NonlinearPoissonData::Field field = NonlinearPoissonData::Field::Sine,
        Boundary condition = Boundary::Dirichlet, bool omitFlux = false,
        bool wrongBoundaryTangent = false)
        : m_mesh(mesh),
          m_order(quadratureOrder),
          m_data(mesh.getDimension(), amplitude, field),
          m_space(std::integral_constant<size_t, K>{}, mesh),
          m_state(m_space),
          m_lift(m_space),
          m_liftBoundary(liftBoundary),
          m_du(m_space),
          m_v(m_space),
          m_problem(m_du, m_v),
          m_ksp(m_problem)
      {
        using namespace Variational;
        m_lift = Zero();
        if (m_liftBoundary)
          m_lift = m_data.getSolution();
        m_state = m_lift;
        const Real gamma = omitCubic ? 0 : 1;
        auto a = Integral(Grad(m_du), Grad(m_v));
        auto c =
          Integral((1 + (wrongTangent ? 1 : 3) * gamma * m_state * m_state) * m_du, m_v);
        auto r = Integral(Grad(m_state), Grad(m_v));
        auto s = Integral(m_state + gamma * m_state * m_state * m_state, m_v);
        auto b = Integral(m_data.getSource(), m_v);
        a.setOrder(m_order);
        c.setOrder(m_order);
        r.setOrder(m_order);
        s.setOrder(m_order);
        b.setOrder(m_order);
        auto boundary = DirichletBC(m_du, Zero());
        if (condition != Boundary::Dirichlet)
          boundary.on(DirichletAttribute);
        if (condition == Boundary::PureNeumann)
          m_freeDOFs = m_space.getSize();
        else
        {
          boundary.assemble();
          const auto& rows =
            std::get<typename Variational::DirichletBCBase<Real>::ValueDOFs>(
              boundary.getDOFs());
          Index begin = 0, end = static_cast<Index>(m_space.getSize());
          if constexpr (requires { mesh.getShard(); })
            m_space.getOwnershipRange(begin, end);
          m_freeDOFs = static_cast<size_t>(end - begin);
          for (const auto& [global, value] : rows)
          {
            (void)value;
            if (global >= begin && global < end)
              --m_freeDOFs;
          }
#ifdef RODIN_USE_MPI
          if constexpr (requires { mesh.getShard(); })
            m_freeDOFs = boost::mpi::all_reduce(
              mesh.getContext().getCommunicator(), m_freeDOFs, std::plus<size_t>());
#endif
        }
        if (condition == Boundary::Dirichlet)
          m_problem = a + c + r + s - b + boundary;
        else
        {
          const BoundaryNormal normal(mesh);
          const size_t dim = mesh.getDimension();
          const RealFunction flux(
            [gradient = m_data.getGradient(), normal, dim](const Geometry::Point& point) {
              const auto derivative = gradient(point);
              if (dim == 1)
                return (point(0) < Real(0.5) ? Real(-1) : Real(1)) * derivative(0);
              const auto outward = normal(point);
              Real value = 0;
              for (size_t j = 0; j < dim; ++j)
                value += derivative(j) * outward(j);
              return value;
            });
          const Real beta = condition == Boundary::Robin ? Real(1) : Real(0);
          auto natural = BoundaryIntegral(
            (omitFlux ? Real(0) : Real(1)) * flux + beta * m_data.getSolution(), m_v);
          natural.setOrder(m_order);
          if (condition == Boundary::PureNeumann)
            m_problem = a + c + r + s - b - natural;
          else if (condition == Boundary::MixedNeumann)
            m_problem = a + c + r + s - b - natural.over(NaturalAttribute) + boundary;
          else
          {
            auto robinResidual = BoundaryIntegral(m_state, m_v);
            auto robinTangent =
              BoundaryIntegral((wrongBoundaryTangent ? Real(0) : Real(1)) * m_du, m_v);
            robinResidual.setOrder(m_order);
            robinTangent.setOrder(m_order);
            m_problem = a + c + r + s - b - natural.over(NaturalAttribute) +
              robinResidual.over(NaturalAttribute) + robinTangent.over(NaturalAttribute) +
              boundary;
          }
        }
        m_problem.assemble();
      }

      /** @brief Global unconstrained dimension, counted by unique DOF ownership. */
      size_t getFreeDOFCount() const
      {
        return m_freeDOFs;
      }

      /** @brief Solve and integrate global manufactured field errors.
       * @note MPI state synchronization, SNES/KSP, residuals, and error norms
       * require every rank of the mesh communicator to participate.
       */
      ErrorNorms solve(Real tolerance = 1e-11)
      {
        return solve(tolerance, 0, [](const auto&, const auto&) {});
      }

      /** @brief Inspect the synchronized converged state with independent norm order. */
      template <class Observer>
      ErrorNorms solve(Real tolerance, size_t normOrder, Observer&& observe)
      {
        Solver::SNES snes(m_ksp);
        snes.setStateUpdate([&](const PETSc::Math::Vector& x) { updateState(x); });
        snes.setType(SNESNEWTONLS).setTolerances(tolerance, tolerance, 1e-14, 20, 1000);
        // SNES uses the KSP handle directly, rather than invoking KSP::solve.
        EXPECT_EQ(KSPSetType(m_ksp.getHandle(), KSPCG), PETSC_SUCCESS);
        EXPECT_EQ(KSPSetTolerances(
                    m_ksp.getHandle(), tolerance * 0.01, tolerance * 0.001, 1e5, 50000),
          PETSC_SUCCESS);
        ::PC pc = nullptr;
        EXPECT_EQ(KSPGetPC(m_ksp.getHandle(), &pc), PETSC_SUCCESS);
        EXPECT_EQ(PCSetType(pc, PCJACOBI), PETSC_SUCCESS);
        PetscReal initialResidual = 0, finalResidual = 0;
        EXPECT_EQ(
          VecNorm(m_problem.getLinearSystem().getVector(), NORM_2, &initialResidual),
          PETSC_SUCCESS);
        if (!m_liftBoundary && m_freeDOFs > 0)
        {
          EXPECT_GT(initialResidual, 0);
        }
        snes.solve();
        EXPECT_TRUE(snes.converged());
        if (!m_liftBoundary && m_freeDOFs > 0)
        {
          EXPECT_GT(snes.getIterationNumber(), 0);
        }
        if (m_freeDOFs == 0)
        {
          EXPECT_EQ(initialResidual, 0);
          EXPECT_EQ(snes.getIterationNumber(), 0);
        }
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
        EXPECT_EQ(KSPGetConvergedReason(m_ksp.getHandle(), &reason), PETSC_SUCCESS);
        if (snes.getIterationNumber() > 0)
        {
          EXPECT_GT(reason, 0);
        }
        // Recompute after synchronizing explicitly, independently of the SNES cache.
        updateState(m_problem.getLinearSystem().getSolution());
        m_problem.assemble(Variational::AssemblyTarget::RHS);
        EXPECT_EQ(
          VecNorm(m_problem.getLinearSystem().getVector(), NORM_2, &finalResidual),
          PETSC_SUCCESS);
        EXPECT_TRUE(std::isfinite(finalResidual));
        EXPECT_LT(finalResidual / std::max(Real(1), initialResidual), 1e-10);
        observe(std::as_const(m_state), std::as_const(m_data));
        return ErrorNorm::compute(m_mesh.get(), m_state, m_data.getSolution(),
          m_data.getGradient(), normOrder == 0 ? m_order : normOrder);
      }

      /** @brief Global residual/Jacobian consistency through SNES callbacks.
       * @note All communicator ranks participate in the perturbed assemblies
       * and vector norms. Directions must vanish on the essential partition.
       */
      Real tangentDefect()
      {
        return tangentDefect(m_data.getSolution());
      }

      template <class Direction>
      Real tangentDefect(const Direction& exactDirection)
      {
        if (m_liftBoundary)
          m_state = Variational::Zero();
        else
          m_state = 0.5 * m_data.getSolution();
        StateType direction(m_space);
        direction = 0.25 * exactDirection;
        Solver::SNES snes(m_ksp);
        snes.setStateUpdate([&](const PETSc::Math::Vector& x) { updateState(x); });
        ::Vec baseline = nullptr, plus = nullptr, minus = nullptr;
        ::Vec action = nullptr, difference = nullptr, residual = nullptr;
        const auto layout = m_problem.getLinearSystem().getSolution();
        EXPECT_EQ(VecDuplicate(layout, &baseline), PETSC_SUCCESS);
        EXPECT_EQ(VecDuplicate(layout, &plus), PETSC_SUCCESS);
        EXPECT_EQ(VecDuplicate(layout, &minus), PETSC_SUCCESS);
        EXPECT_EQ(VecDuplicate(layout, &action), PETSC_SUCCESS);
        EXPECT_EQ(VecDuplicate(layout, &difference), PETSC_SUCCESS);
        EXPECT_EQ(VecDuplicate(layout, &residual), PETSC_SUCCESS);
        EXPECT_EQ(VecCopy(m_state.getData(), baseline), PETSC_SUCCESS);
        constexpr Real epsilon = 1e-5;
        EXPECT_EQ(VecCopy(baseline, plus), PETSC_SUCCESS);
        EXPECT_EQ(VecCopy(baseline, minus), PETSC_SUCCESS);
        EXPECT_EQ(VecAXPY(plus, epsilon, direction.getData()), PETSC_SUCCESS);
        EXPECT_EQ(VecAXPY(minus, -epsilon, direction.getData()), PETSC_SUCCESS);
        auto matrix = m_problem.getLinearSystem().getOperator();
        EXPECT_EQ(
          SNESComputeJacobian(snes.getHandle(), baseline, matrix, matrix), PETSC_SUCCESS);
        EXPECT_EQ(MatMult(matrix, direction.getData(), action), PETSC_SUCCESS);
        EXPECT_EQ(SNESComputeFunction(snes.getHandle(), plus, difference), PETSC_SUCCESS);
        EXPECT_EQ(SNESComputeFunction(snes.getHandle(), minus, residual), PETSC_SUCCESS);
        EXPECT_EQ(VecAXPY(difference, -1, residual), PETSC_SUCCESS);
        EXPECT_EQ(VecScale(difference, 1 / (2 * epsilon)), PETSC_SUCCESS);
        PetscReal denominator = 0, defect = 0;
        EXPECT_EQ(VecNorm(difference, NORM_2, &denominator), PETSC_SUCCESS);
        EXPECT_TRUE(std::isfinite(denominator));
        EXPECT_GT(denominator, 0);
        EXPECT_EQ(VecAXPY(action, -1, difference), PETSC_SUCCESS);
        EXPECT_EQ(VecNorm(action, NORM_2, &defect), PETSC_SUCCESS);
        EXPECT_EQ(VecDestroy(&residual), PETSC_SUCCESS);
        EXPECT_EQ(VecDestroy(&difference), PETSC_SUCCESS);
        EXPECT_EQ(VecDestroy(&action), PETSC_SUCCESS);
        EXPECT_EQ(VecDestroy(&minus), PETSC_SUCCESS);
        EXPECT_EQ(VecDestroy(&plus), PETSC_SUCCESS);
        EXPECT_EQ(VecDestroy(&baseline), PETSC_SUCCESS);
        return defect / denominator;
      }

    private:
      void updateState(const PETSc::Math::Vector& x)
      {
        m_state.setData(x);
        if (m_liftBoundary)
          m_state.axpy(1, m_lift);
      }
      std::reference_wrapper<const MeshType> m_mesh;
      size_t m_order;
      NonlinearPoissonData m_data;
      SpaceType m_space;
      StateType m_state;
      StateType m_lift;
      bool m_liftBoundary;
      size_t m_freeDOFs = 0;
      TrialType m_du;
      TestType m_v;
      ProblemType m_problem;
      Solver::KSP m_ksp;
  };
}

#endif
