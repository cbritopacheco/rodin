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
#include "../convergence/CurvedGeometry.h"
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
   * nonlinear iterations, or linear solves in this scope. Curved Q2 cases use
   * the existing volume-preserving quadratic shear in dimensions two and
   * three. Its first coordinate is unchanged, so the same affine state and
   * analytic actions apply without a geometry-dependent reference solve.
   */
  template <size_t K>
  class NonlinearPoissonAssembly
  {
    public:
      static void run(benchmark::State& timing, Geometry::Polytope::Type geometry,
        bool residual, bool curved)
      {
        using namespace Variational;
        const size_t n = timing.range(0), order = timing.range(1);
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(n);
        Optional<Convergence::CurvedGeometry<decltype(mesh)>> map;
        if (curved)
        {
          assert(mesh.getDimension() >= 2);
          map.emplace(mesh);
          map->template install<2>();
        }
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
        timing.counters["geometry_degree"] = curved ? 2 : 1;
        timing.counters["quadrature_order"] = order;
        timing.counters["quadrature_points_per_cell"] =
          QF::PolytopeQuadratureFormula::get(order, geometry).getSize();
      }

      /**
       * @brief Measures a source load that requires physical coordinates.
       *
       * The timed form is @f$L_s(v)=s\int_\Omega(1+x_0+x_{d-1}^2)v\,dx@f$.
       * An independent pointwise map, Jacobian and basis loop checks every
       * assembled coefficient. Constant-test actions have closed-form moments
       * on the unit box and its prescribed quadratic image. Unlike the affine
       * state benchmarks, this source study also supports the curved segment.
       *
       * @param timing Benchmark iteration state and workload sizes.
       * @param geometry Reference-cell family.
       * @param curved Whether to install exact quadratic geometry.
       */
      static void runSource(
        benchmark::State& timing, Geometry::Polytope::Type geometry, bool curved)
      {
        using namespace Variational;
        const size_t n = timing.range(0), order = timing.range(1);
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(n);
        const size_t dimension = mesh.getDimension();
        constexpr Real Amplitude = Real(1) / 10;
        Optional<Convergence::CurvedGeometry<decltype(mesh)>> map;
        if (curved)
        {
          map.emplace(
            mesh, Convergence::CurvedGeometry<decltype(mesh)>::Map::Quadratic, Amplitude);
          map->template install<2>();
        }
        H1 space(std::integral_constant<size_t, K>{}, mesh);
        GridFunction one(space);
        one = RealFunction(1);
        TestFunction v(space);
        Real scale = 1;
        const RealFunction source([&](const Geometry::Point& point) {
          const auto& x = point.getPhysicalCoordinates();
          return scale * (1 + x(0) + x(dimension - 1) * x(dimension - 1));
        });
        auto integral = Integral(source, v);
        integral.setOrder(order);
        LinearForm load(v);
        load = integral;

        Math::Vector<Real> reference = Math::Vector<Real>::Zero(space.getSize());
        const auto& formula = QF::PolytopeQuadratureFormula::get(order, geometry);
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const auto& element = space.getFiniteElement(dimension, cell->getIndex());
          const auto dofs = space.getDOFs(dimension, cell->getIndex());
          for (size_t qp = 0; qp < formula.getSize(); ++qp)
          {
            const Geometry::Point point(*cell, formula.getPoint(qp));
            const auto& x = point.getPhysicalCoordinates();
            const Real weight = formula.getWeight(qp) * point.getDistortion();
            const Real value = 1 + x(0) + x(dimension - 1) * x(dimension - 1);
            for (size_t local = 0; local < element.getCount(); ++local)
            {
              reference(dofs(local)) +=
                weight * value * element.getBasis(local)(formula.getPoint(qp));
            }
          }
        }
        const Real amplitude = curved ? Amplitude : 0;
        const Real length = 1 + amplitude;
        const Real moment = dimension == 1
          ? length + length * length / 2 + length * length * length / 3
          : Real(11) / 6 + amplitude / 3 + amplitude * amplitude / 5;
        const Real omittedMoment =
          dimension == 1 ? length + length * length / 2 : Real(3) / 2;
        constexpr Real ActionTolerance = 1e-9;
        constexpr Real VectorTolerance = 1e-11;
        constexpr Real MissingQuadraticMinimumDefect = 0.1;
        const auto valid = [&] {
          const auto& values = load.getVector();
          const Real action = one.getData().dot(values);
          return values.allFinite() && std::isfinite(action) &&
            std::abs(action / (scale * moment) - 1) < ActionTolerance &&
            (values - scale * reference).norm() <
            VectorTolerance * std::max(Real(1), (scale * reference).norm()) &&
            std::abs(omittedMoment / moment - 1) > MissingQuadraticMinimumDefect;
        };
        load.assemble();
        if (!valid())
        {
          nonlinearCheckFailed = true;
          timing.SkipWithError("Physical source vector or analytic moment failed");
          return;
        }
        for (auto _ : timing)
        {
          load.assemble();
          benchmark::DoNotOptimize(load.getVector().data());
          benchmark::ClobberMemory();
        }
        bool ok = valid();
        scale = 2;
        load.assemble();
        ok = ok && valid();
        if (!ok)
        {
          nonlinearCheckFailed = true;
          timing.SkipWithError("Physical source replacement or state update failed");
        }
        timing.counters["cells"] = mesh.getPolytopeCount(dimension);
        timing.counters["dofs"] = space.getSize();
        timing.counters["degree"] = K;
        timing.counters["geometry_degree"] = curved ? 2 : 1;
        timing.counters["quadrature_order"] = order;
        timing.counters["quadrature_points_per_cell"] = formula.getSize();
      }

      static void registerCases()
      {
        using G = Geometry::Polytope::Type;
        for (auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
               G::Hexahedron, G::Wedge, G::Pyramid})
        {
          for (bool curved : {false, true})
          {
            const std::string name = std::string("NonlinearPoisson/") +
              (curved ? "CurvedQ2/" : "") + "Source/P" + std::to_string(K) + "/" +
              std::string(Convergence::UniformGrid::getGeometryName(geometry));
            benchmark::RegisterBenchmark(name.c_str(),
              [geometry, curved](auto& state) { runSource(state, geometry, curved); })
              ->ArgsProduct({{3, 5, 9}, {8, 16}})
              ->UseRealTime();
          }
          for (bool residual : {false, true})
          {
            for (bool curved : {false, true})
            {
              // In one dimension the shear changes the volume and affine pullback.
              if (curved && geometry == G::Segment)
                continue;
              const std::string name = std::string("NonlinearPoisson/") +
                (curved ? "CurvedQ2/" : "") + (residual ? "Residual" : "Tangent") + "/P" +
                std::to_string(K) + "/" +
                std::string(Convergence::UniformGrid::getGeometryName(geometry));
              benchmark::RegisterBenchmark(name.c_str(),
                [geometry, residual, curved](
                  auto& state) { run(state, geometry, residual, curved); })
                ->ArgsProduct({{3, 5, 9}, {8, 16}})
                ->UseRealTime();
            }
          }
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
