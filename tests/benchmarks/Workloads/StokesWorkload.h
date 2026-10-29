/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_WORKLOADS_STOKESWORKLOAD_H
#define RODIN_WORKLOADS_STOKESWORKLOAD_H

#include <benchmark/benchmark.h>
#include <sstream>
#include "Rodin/Assembly.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Solver/GMRES.h"
#include "StokesBlockPreconditioner.h"
#include "convergence/Stokes.h"

namespace Rodin::Tests::Benchmarks
{
  /**
   * @brief Shared native Taylor--Hood assembly and fixed-system solve workload.
   *
   * ## Mathematical contract
   *
   * On the unit box, @f$u=x_1e_0@f$, @f$p=x_0-1/2@f$ and @f$f=e_0@f$
   * solve @f$-\Delta u+\nabla p=f@f$ with zero divergence and pressure mean.
   * Velocity traces are prescribed. A global P0g multiplier enforces the
   * pressure mean, as in Stokes convergence certification. Trial/test block
   * offsets are obtained from Problem; no interleaved DOF layout is assumed.
   *
   * ## Measurement boundary
   *
   * Assembly measures public repeated Problem::assemble, including constraints
   * and global completion. Solve measures public SparseLU::solve(system),
   * including fresh factorization and one residual correction, but no assembly
   * or field scatter. It does not claim to isolate triangular substitution.
   * Block GMRES instead measures preconditioner setup and the restarted solve
   * from zero, without a correction solve. Pressure mass assembly is untimed.
   * Assembly qualification does not construct or run a solver.
   * Setup and all numerical qualifications are outside timing. Both variants
   * check the exact affine patch residual, a wrong-field control and a
   * missing-pressure-gradient control. The solver variant additionally checks
   * solver success, coefficient reproduction and the true solved residual.
   *
   * @tparam K Velocity degree; pressure has degree K-1.
   * @tparam Weighted Whether explicit unit viscosity selects the weighted kernel.
   * @tparam BlockGMRES Whether GMRES with the fixed block preconditioner replaces SparseLU.
   */
  template <size_t K, bool Weighted = false, bool BlockGMRES = false>
  class StokesWorkload
  {
    public:
      /// @brief Common order used by the certified mixed polynomial workload.
      static constexpr size_t QuadratureOrder = 16;

      /**
       * @brief Defines the certified mixed form or a missing-pressure control.
       * @param problem Destination problem.
       * @param u Velocity trial function.
       * @param p Pressure trial function.
       * @param lambda Global mean multiplier trial function.
       * @param v Velocity test function.
       * @param q Pressure test function.
       * @param mu Global mean test function.
       * @param data Analytic fields and forcing.
       * @param omitPressure Whether the momentum pressure gradient is omitted.
       */
      template <class ProblemType, class U, class P, class L, class V, class Q, class M>
      static void configure(ProblemType& problem, U& u, P& p, L& lambda, V& v, Q& q,
        M& mu, const Convergence::StokesData& data, bool omitPressure = false)
      {
        using namespace Variational;
        auto diffusion = [&]() {
          if constexpr (Weighted)
            return Integral(Real(1) * Jacobian(u), Jacobian(v));
          else
            return Integral(Jacobian(u), Jacobian(v));
        }();
        auto momentum = Integral(p, Div(v));
        auto divergence = Integral(Div(u), q);
        auto gaugePressure = Integral(lambda, q);
        auto gaugeMean = Integral(p, mu);
        auto body = Integral(data.getForcing(), v);
        diffusion.setOrder(QuadratureOrder);
        momentum.setOrder(QuadratureOrder);
        divergence.setOrder(QuadratureOrder);
        gaugePressure.setOrder(QuadratureOrder);
        gaugeMean.setOrder(QuadratureOrder);
        body.setOrder(QuadratureOrder);
        if (omitPressure)
        {
          problem = diffusion + divergence + gaugePressure + gaugeMean - body +
            DirichletBC(u, data.getVelocity());
        }
        else
        {
          problem = diffusion - momentum + divergence + gaugePressure + gaugeMean - body +
            DirichletBC(u, data.getVelocity());
        }
      }

      /**
       * @brief Measures one public operation on a fixed mixed discretization.
       * @tparam Solve Whether fixed-system solving replaces assembly in timing.
       * @param state Points per axis and repetitions.
       * @param geometry Unit-box cell family.
       * @param failed Process-level qualification failure flag.
       */
      template <bool Solve>
      static void run(
        benchmark::State& state, Geometry::Polytope::Type geometry, bool& failed)
      {
        using namespace Variational;
        static_assert(K >= 2);
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
        const size_t dim = mesh.getDimension();
        const Convergence::StokesData data(dim, Convergence::StokesData::Field::Affine);
        H1 velocitySpace(std::integral_constant<size_t, K>{}, mesh, dim);
        H1 pressureSpace(std::integral_constant<size_t, K - 1>{}, mesh);
        P0g meanSpace(mesh);
        TrialFunction u(velocitySpace);
        TrialFunction p(pressureSpace);
        TrialFunction lambda(meanSpace);
        TestFunction v(velocitySpace);
        TestFunction q(pressureSpace);
        TestFunction mu(meanSpace);
        Problem problem(u, p, lambda, v, q, mu), omitted(u, p, lambda, v, q, mu);
        configure(problem, u, p, lambda, v, q, mu, data);
        configure(omitted, u, p, lambda, v, q, mu, data, true);
        u.getSolution() = data.getVelocity();
        p.getSolution() = data.getPressure();
        lambda.getSolution() = Real(0);
        problem.assemble();
        omitted.assemble();
        auto& system = problem.getLinearSystem();
        const Math::Vector<Real> exact = system.getSolution();
        const auto& baseline = omitted.getLinearSystem();
        const Real rhsScale = std::max(Real(1), system.getVector().norm());
        const Real exactScale = std::max(Real(1), exact.norm());
        const auto residual = [&](const auto& x) {
          return (system.getOperator() * x - system.getVector()).norm() / rhsScale;
        };
        const Real patchResidual = residual(exact);
        const Real omittedResidual =
          (baseline.getOperator() * exact - baseline.getVector()).norm() / rhsScale;
        const Math::Vector<Real> wrong = exact.array() + Real(1);
        const Real wrongResidual = residual(wrong);
        using System = std::remove_reference_t<decltype(system)>;
        using SolverType = std::conditional_t<BlockGMRES,
          Solver::GMRES<System, StokesBlockPreconditioner>, Solver::SparseLU<System>>;
        Optional<SolverType> solver;
        if constexpr (Solve)
        {
          solver.emplace(problem);
          if constexpr (BlockGMRES)
          {
            BilinearForm pressureMass(p, q);
            auto integrand = Integral(p, q);
            integrand.setOrder(QuadratureOrder);
            pressureMass = integrand;
            pressureMass.assemble();
            solver->getPreconditioner().configure(
              std::move(pressureMass.getOperator()), problem.getTrialOffsets()[1]);
            solver->setRestart(60).setMaxIterations(1200).setTolerance(1e-13);
            system.getSolution().setZero();
          }
          else
            solver->setRefinementSteps(1);
          solver->solve(system);
        }
        const auto coefficientError = [&]() {
          return (system.getSolution() - exact).norm() / exactScale;
        };
        const auto valid = [&]() {
          if constexpr (Solve)
          {
            if (!solver->success() || !std::isfinite(coefficientError()) ||
              coefficientError() >= 1e-9 ||
              !std::isfinite(residual(system.getSolution())) ||
              residual(system.getSolution()) >= 1e-11)
              return false;
          }
          return std::isfinite(patchResidual) && patchResidual < 1e-11 &&
            std::isfinite(omittedResidual) && omittedResidual > 1e-4 &&
            std::isfinite(wrongResidual) && wrongResidual > 1e-4;
        };
        if (!valid())
        {
          failed = true;
          std::ostringstream message;
          message << "Stokes qualification failed: patch=" << patchResidual
                  << ", missing-pressure=" << omittedResidual
                  << ", wrong-field=" << wrongResidual;
          if constexpr (Solve)
          {
            message << ", solver-success=" << solver->success()
                    << ", coefficient=" << coefficientError()
                    << ", solved-residual=" << residual(system.getSolution());
          }
          state.SkipWithError(message.str().c_str());
          return;
        }
        bool repeatedValid = true;
        for (auto _ : state)
        {
          if constexpr (Solve)
          {
            if constexpr (BlockGMRES)
              system.getSolution().setZero();
            solver->solve(system);
          }
          else
            problem.assemble();
          benchmark::DoNotOptimize(system.getSolution().data());
          benchmark::ClobberMemory();
          state.PauseTiming();
          if constexpr (Solve)
            repeatedValid = repeatedValid && valid();
          else
            repeatedValid = repeatedValid && residual(exact) < 1e-11;
          state.ResumeTiming();
        }
        if (!repeatedValid)
        {
          failed = true;
          state.SkipWithError("Repeated Stokes operation changed numerical results");
        }
        state.counters["dimension"] = dim;
        state.counters["velocity_degree"] = K;
        state.counters["pressure_degree"] = K - 1;
        state.counters["cells"] = mesh.getPolytopeCount(dim);
        state.counters["velocity_dofs"] = velocitySpace.getSize();
        state.counters["pressure_dofs"] = pressureSpace.getSize();
        state.counters["mean_dofs"] = meanSpace.getSize();
        state.counters["explicit_unit_viscosity"] = Weighted;
        state.counters["rows"] = system.getOperator().rows();
        state.counters["nnz"] = system.getOperator().nonZeros();
        state.counters["quadrature_order"] = QuadratureOrder;
        state.counters["quadrature_points_per_cell"] =
          QF::QuadratureFormula::get(QuadratureOrder, geometry).getSize();
        state.counters["relative_patch_residual"] = patchResidual;
        state.counters["relative_missing_pressure_residual"] = omittedResidual;
        if constexpr (Solve)
        {
          if constexpr (BlockGMRES)
          {
            state.counters["iterations"] = solver->getIterationNumber();
            state.counters["estimated_relative_error"] = solver->getError();
            state.counters["restart"] = 60;
          }
          state.counters["relative_coefficient_error"] = coefficientError();
          state.counters["relative_residual"] = residual(system.getSolution());
        }
      }

      /**
       * @brief Registers six applicable cell families and three mesh sizes.
       * @tparam Solve Whether to register solving rather than assembly.
       * @param failed Process-level qualification failure flag.
       */
      template <bool Solve>
      static void registerCases(bool& failed)
      {
        using Geometry::Polytope;
        for (const auto geometry :
          {Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
            Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
            Polytope::Type::Hexahedron, Polytope::Type::Wedge})
        {
          const std::string name = std::string("Stokes/") +
            (Solve ? (BlockGMRES ? "GMRESBlockSchur/" : "SparseLU/") : "Assembly/") +
            (Weighted ? "WeightedViscosity/" : "") + "H1/P" + std::to_string(K) + "P" +
            std::to_string(K - 1) + "/Real/" +
            std::string(Convergence::UniformGrid::getGeometryName(geometry));
          benchmark::RegisterBenchmark(name.c_str(),
            [geometry, &failed](auto& state) { run<Solve>(state, geometry, failed); })
            ->Arg(3)
            ->Arg(5)
            ->Arg(9)
            ->UseRealTime();
        }
      }
  };
}
#endif
