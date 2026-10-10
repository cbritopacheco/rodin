/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Prescribed-state semilinear residual and tangent assembly.
 */
#include <benchmark/benchmark.h>

#include "Rodin/Variational/H1.h"
#include "../convergence/Convergence.h"
#include "../QuadratureReference.h"
#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool nonlinearCheckFailed = false;

  /**
   * @brief Measures the volume assembly used by semilinear Poisson Newton steps.
   *
   * The benchmark separates residual and tangent costs at a fixed discrete state.
   *
   * ## Architecture
   *
   * On the unit box, a public DOF interpolation prescribes @f$q=s(1+x_0)@f$.
   * Separate public forms assemble
   * @f$R(q;v)=\int\nabla q\cdot\nabla v+(q+q^3)v@f$ and
   * @f$J(q;w,v)=\int\nabla w\cdot\nabla v+(1+3q^2)wv@f$.
   * Untimed gates compare the tangent with the original quadrature loop,
   * check analytic actions, and differentiate the assembled residual. Warm
   * replacement assembly is timed; a post-timing state change checks cache
   * invalidation. There are no constraints, boundary integrals, source terms,
   * nonlinear iterations, or linear solves in this scope.
   */
  template <size_t K>
  class NonlinearPoissonAssembly
  {
    public:
      static void run(
        benchmark::State& timing, Geometry::Polytope::Type geometry, bool residual)
      {
        using namespace Variational;
        const size_t n = timing.range(0), order = timing.range(1);
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(n);
        H1 space(std::integral_constant<size_t, K>{}, mesh);
        GridFunction q(space);
        TrialFunction u(space);
        TestFunction v(space);
        const RealFunction affine(
          [](const Geometry::Point& point) { return 1 + point(0); });
        q = affine;
        u.getSolution() = affine;
        BilinearForm tangent(u, v), reference(u, v);
        LinearForm load(v);
        auto diffusion = Integral(Grad(u), Grad(v));
        auto reaction = Integral((1 + 3 * q * q) * u, v);
        auto gradientResidual = Integral(Grad(q), Grad(v));
        auto reactionResidual = Integral(q + q * q * q, v);
        diffusion.setOrder(order);
        reaction.setOrder(order);
        gradientResidual.setOrder(order);
        reactionResidual.setOrder(order);
        tangent = diffusion + reaction;
        reference = ReferenceIntegral(diffusion) + ReferenceIntegral(reaction);
        load = gradientResidual + reactionResidual;

        // Scale-relative numerical oracles, not runtime acceptance thresholds.
        constexpr Real ActionTolerance = 1e-9;
        constexpr Real OperatorTolerance = 1e-11;
        constexpr Real DifferenceStep = 1e-5;
        constexpr Real DerivativeTolerance = 1e-7;
        constexpr Real MissingCubicMinimumDefect = 0.1;
        constexpr Real AffineSquareIntegral = Real(7) / 3;
        constexpr Real AffineFourthIntegral = Real(31) / 5;
        const auto valid = [&](Real scale) {
          reference.assemble();
          const auto& coefficients = q.getData();
          const Real square = scale * scale, fourth = square * square;
          const Real residualAction =
            square * (1 + AffineSquareIntegral) + fourth * AffineFourthIntegral;
          const Real tangentAction =
            square * (1 + AffineSquareIntegral) + 3 * fourth * AffineFourthIntegral;
          const Real actualResidual = coefficients.dot(load.getVector());
          const Real actualTangent =
            coefficients.dot(tangent.getOperator() * coefficients);
          const auto& baseline = reference.getOperator();
          return std::isfinite(actualResidual) && std::isfinite(actualTangent) &&
            std::abs(actualResidual / residualAction - 1) < ActionTolerance &&
            std::abs(actualTangent / tangentAction - 1) < ActionTolerance &&
            (tangent.getOperator() - baseline).norm() <
            OperatorTolerance * std::max(Real(1), baseline.norm()) &&
            std::abs((square * (1 + AffineSquareIntegral)) / residualAction - 1) >
            MissingCubicMinimumDefect;
        };
        tangent.assemble();
        load.assemble();
        bool ok = valid(1);

        const auto direction = u.getSolution().getData();
        q = (1 + DifferenceStep) * affine;
        load.assemble();
        const auto plus = load.getVector();
        q = (1 - DifferenceStep) * affine;
        load.assemble();
        const auto difference = ((plus - load.getVector()) / (2 * DifferenceStep)).eval();
        const auto action = (tangent.getOperator() * direction).eval();
        ok = ok && difference.allFinite() &&
          (difference - action).norm() <
            DerivativeTolerance * std::max(Real(1), action.norm());
        q = affine;
        tangent.assemble();
        load.assemble();
        ok = ok && valid(1);
        if (!ok)
        {
          nonlinearCheckFailed = true;
          timing.SkipWithError(
            "Nonlinear action, reference operator, or tangent oracle failed");
          return;
        }
        for (auto _ : timing)
        {
          if (residual)
          {
            load.assemble();
            benchmark::DoNotOptimize(load.getVector().data());
          }
          else
          {
            tangent.assemble();
            benchmark::DoNotOptimize(tangent.getOperator().nonZeros());
          }
          benchmark::ClobberMemory();
        }
        ok = valid(1);
        q = 2 * affine;
        tangent.assemble();
        load.assemble();
        ok = ok && valid(2);
        if (!ok)
        {
          nonlinearCheckFailed = true;
          timing.SkipWithError(
            "Repeated assembly or nonlinear state update changed numerics");
        }
        timing.counters["cells"] = mesh.getPolytopeCount(mesh.getDimension());
        timing.counters["dofs"] = space.getSize();
        timing.counters["nnz"] = tangent.getOperator().nonZeros();
        timing.counters["degree"] = K;
        timing.counters["quadrature_order"] = order;
        timing.counters["quadrature_points_per_cell"] =
          QF::PolytopeQuadratureFormula::get(order, geometry).getSize();
      }

      static void registerCases()
      {
        using G = Geometry::Polytope::Type;
        for (auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
               G::Hexahedron, G::Wedge, G::Pyramid})
          for (bool residual : {false, true})
          {
            const std::string name = std::string("NonlinearPoisson/") +
              (residual ? "Residual" : "Tangent") + "/P" + std::to_string(K) + "/" +
              std::string(Convergence::UniformGrid::getGeometryName(geometry));
            benchmark::RegisterBenchmark(name.c_str(),
              [geometry, residual](auto& state) { run(state, geometry, residual); })
              ->ArgsProduct({{3, 5, 9}, {8, 16}})
              ->UseRealTime();
          }
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  NonlinearPoissonAssembly<1>::registerCases();
  NonlinearPoissonAssembly<2>::registerCases();
  NonlinearPoissonAssembly<3>::registerCases();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
#ifdef RODIN_USE_OPENMP
  benchmark::AddCustomContext("assembly_backend", "Eigen/OpenMP");
  benchmark::AddCustomContext(
    "openmp_max_threads", std::to_string(omp_get_max_threads()));
#else
  benchmark::AddCustomContext("assembly_backend", "Eigen/Sequential");
#endif
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t matched = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return nonlinearCheckFailed || matched == 0 ? 1 : 0;
}
