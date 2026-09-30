/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Warmed real P1/P2 physics operator assembly on all geometries. */

#include <benchmark/benchmark.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational/H1.h"
#include "../convergence/Convergence.h"

#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  bool assemblyCheckFailed = false;

  /**
   * @brief Measures repeated assembly of a fixed physical bilinear form.
   *
   * Architecture: mesh, space, and form construction precede timing. An
   * interpolated affine field supplies an independent analytic energy oracle.
   * Two untimed assemblies verify replacement rather than accumulation. The
   * timed loop includes the complete public assembly operation; the final
   * operator is checked outside timing against the initial operator and energy.
   * No solve, constraint elimination, or mesh construction is timed.
   */
  template <size_t Degree, bool Elasticity>
  class PhysicsAssembly
  {
    public:
      static void run(
        benchmark::State& state, Polytope::Type geometry, bool conductivity = false)
      {
        const size_t n = state.range(0);
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(n);
        const size_t dim = mesh.getDimension();
        using Range = std::conditional_t<Elasticity, Math::SpatialVector<Real>, Real>;
        auto space = [&]() {
          if constexpr (Elasticity)
            return H1<Degree, Range, LocalMesh>(
              std::integral_constant<size_t, Degree>{}, mesh, dim);
          else
            return H1<Degree, Range, LocalMesh>(
              std::integral_constant<size_t, Degree>{}, mesh);
        }();
        TrialFunction u(space);
        TestFunction v(space);
        BilinearForm form(u, v);
        Real expected;
        if constexpr (Elasticity)
        {
          auto volumetric = Integral(1.5 * Div(u), Div(v));
          auto shear = Integral(
            0.5 * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()));
          volumetric.setOrder(6);
          shear.setOrder(6);
          form = volumetric + shear;
          u.getSolution() = VectorFunction(dim, [dim](const Point& p) {
            Real s = 0;
            for (size_t j = 0; j < dim; ++j)
              s += p(j);
            Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
            for (size_t i = 0; i < dim; ++i)
              value(i) = Real(i + 1) * s;
            return value;
          });
          const Real c = Real(dim * (dim + 1)) / 2;
          const Real squares = Real(dim * (dim + 1) * (2 * dim + 1)) / 6;
          expected = 1.5 * c * c + 0.5 * (Real(dim) * squares + c * c);
        }
        else
        {
          const RealFunction gamma([dim, conductivity](const Point& p) {
            Real value = 1;
            if (conductivity)
              for (size_t j = 0; j < dim; ++j)
                value += p(j);
            return value;
          });
          auto diffusion = Integral(gamma * Grad(u), Grad(v));
          diffusion.setOrder(6);
          form = diffusion;
          u.getSolution() = RealFunction([dim](const Point& p) {
            Real value = 0;
            for (size_t j = 0; j < dim; ++j)
              value += p(j);
            return value;
          });
          expected = Real(dim) * (conductivity ? 1 + Real(dim) / 2 : 1);
        }

        form.assemble();
        const auto baseline = form.getOperator();
        const auto& coefficients = u.getSolution().getData();
        const auto valid = [&]() {
          const auto& matrix = form.getOperator();
          const Real energy = coefficients.dot(matrix * coefficients);
          return std::isfinite(energy) && std::abs(energy / expected - 1) < 1e-9 &&
            (matrix - baseline).norm() < 1e-12 * std::max(Real(1), baseline.norm());
        };
        form.assemble();
        if (!valid())
        {
          assemblyCheckFailed = true;
          state.SkipWithError("Assembly failed analytic-energy or replacement check");
          return;
        }

        for (auto _ : state)
        {
          form.assemble();
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
          benchmark::ClobberMemory();
        }
        if (!valid())
        {
          assemblyCheckFailed = true;
          state.SkipWithError("Timed assembly changed operator or analytic energy");
        }
        const size_t cells = mesh.getPolytopeCount(dim);
        state.counters["cells"] = cells;
        state.counters["dofs"] = space.getSize();
        state.counters["nnz"] = form.getOperator().nonZeros();
        state.counters["degree"] = Degree;
        state.counters["quadrature_order"] = 6;
        state.counters["components"] = Elasticity ? dim : 1;
        state.counters["cells_per_second"] =
          benchmark::Counter(cells, benchmark::Counter::kIsIterationInvariantRate);
      }
  };

  /// @brief Registers three mesh sizes for each physical context and geometry.
  template <size_t Degree>
  void registerPhysics()
  {
    for (auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      const std::string suffix = "/P" + std::to_string(Degree) + "/" +
        std::string(Convergence::UniformGrid::getGeometryName(geometry));
      for (bool conductivity : {false, true})
      {
        const std::string name =
          std::string(conductivity ? "Conductivity" : "Poisson") + suffix;
        benchmark::RegisterBenchmark(name.c_str(),
          [geometry, conductivity](auto& state) {
            PhysicsAssembly<Degree, false>::run(state, geometry, conductivity);
          })
          ->Arg(3)
          ->Arg(5)
          ->Arg(9)
          ->UseRealTime();
      }
      const std::string name = "LinearElasticity" + suffix;
      benchmark::RegisterBenchmark(name.c_str(),
        [geometry](auto& state) { PhysicsAssembly<Degree, true>::run(state, geometry); })
        ->Arg(3)
        ->Arg(5)
        ->Arg(9)
        ->UseRealTime();
    }
  }
}

int main(int argc, char** argv)
{
  Rodin::Tests::Benchmarks::registerPhysics<1>();
  Rodin::Tests::Benchmarks::registerPhysics<2>();
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
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return Rodin::Tests::Benchmarks::assemblyCheckFailed ? 1 : 0;
}
