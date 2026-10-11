/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Local profiling of the certified finite-mesh pressure-spectrum oracle.
 */
#include <benchmark/benchmark.h>
#include <array>
#include <string>
#include "Rodin/Assembly.h"
#include "../convergence/Convergence.h"
#include "../convergence/CurvedGeometry.h"
#include "../convergence/MixedStability.h"

namespace Rodin::Tests::Benchmarks
{
  using namespace Geometry;
  using namespace Variational;
  using namespace Convergence;

  bool pressureOracleFailed = false;

  /**
   * @brief Measures the complete oracle on fixed Taylor--Hood operators.
   *
   * The velocity energy, divergence and pressure mass matrices define the
   * zero-mean generalized Schur spectrum. Both the Schur eigensolver and the
   * independently whitened-divergence singular-value route are retained.
   *
   * ## Architecture
   *
   * Mesh construction, FE assembly, logical boundary constraints and pressure
   * constant interpolation precede timing. Every timed invocation constructs
   * a fresh result and includes reduction, factorization, blocked solves, QR,
   * spectra and residual diagnostics. Untimed warmup/repeat checks guard the
   * measurement; they do not count as additional convergence certification.
   * The global oracle is rank-local once supplied with complete matrices.
   * This benchmark introduces no MPI operations and measures no sharding.
   *
   * @tparam K Velocity degree; pressure degree is one less.
   */
  template <size_t K>
  class PressureOracle
  {
    public:
      /**
       * @brief Assembles fixed inputs and profiles repeated oracle evaluation.
       * @param state Grid points, quadrature order and benchmark iterations.
       * @param geometry Reference-cell family.
       * @param curved Whether to install the exact quadratic map.
       */
      static void run(benchmark::State& state, Polytope::Type geometry, bool curved)
      {
        auto mesh = UniformGrid(geometry).makeMesh(state.range(0));
        if (curved)
        {
          CurvedGeometry mapping(mesh);
          mapping.template install<2>();
        }
        H1 velocity(std::integral_constant<size_t, K>{}, mesh, mesh.getDimension());
        H1 pressure(std::integral_constant<size_t, K - 1>{}, mesh);
        TrialFunction u(velocity);
        TrialFunction p(pressure);
        TestFunction v(velocity);
        TestFunction q(pressure);
        auto diffusion = Integral(Jacobian(u), Jacobian(v));
        auto divergence = Integral(Div(u), q);
        auto pressureMass = Integral(p, q);
        const size_t order = state.range(1);
        diffusion.setOrder(order);
        divergence.setOrder(order);
        pressureMass.setOrder(order);
        BilinearForm a(u, v);
        BilinearForm b(u, q);
        BilinearForm mass(p, q);
        a = diffusion;
        b = divergence;
        mass = pressureMass;
        a.assemble();
        b.assemble();
        mass.assemble();
        auto boundary = DirichletBC(u, Zero());
        boundary.assemble();
        GridFunction constant(pressure);
        constant = RealFunction(Real(1));
        const auto& constraints = std::get<IndexMap<Real>>(boundary.getDOFs());
        MixedStability::Result baseline;
        MixedStability::compute(baseline, a.getOperator(), b.getOperator(),
          mass.getOperator(), constraints, constant.getData());
        if (!::testing::Test::HasFailure())
          MixedStability::expectConsistent(baseline);
        if (::testing::Test::HasFailure())
        {
          pressureOracleFailed = true;
          state.SkipWithError("Pressure-oracle warmup failed");
          return;
        }
        for (auto _ : state)
        {
          MixedStability::Result measured;
          MixedStability::compute(measured, a.getOperator(), b.getOperator(),
            mass.getOperator(), constraints, constant.getData());
          benchmark::DoNotOptimize(measured.eigenvalues.data());
          benchmark::ClobberMemory();
          if (::testing::Test::HasFailure())
          {
            pressureOracleFailed = true;
            state.SkipWithError("Timed pressure oracle failed");
            break;
          }
        }
        MixedStability::Result final;
        MixedStability::compute(final, a.getOperator(), b.getOperator(),
          mass.getOperator(), constraints, constant.getData());
        if (!::testing::Test::HasFailure())
          MixedStability::expectConsistent(final);
        if (::testing::Test::HasFailure() ||
          final.eigenvalues.size() != baseline.eigenvalues.size() ||
          !(final.eigenvalues.array() == baseline.eigenvalues.array()).all())
        {
          pressureOracleFailed = true;
          state.SkipWithError("Repeated pressure oracle changed its spectrum");
        }
        state.counters["dimension"] = mesh.getDimension();
        state.counters["cells"] = mesh.getPolytopeCount(mesh.getDimension());
        state.counters["velocity_dofs"] = velocity.getSize();
        state.counters["pressure_dofs"] = pressure.getSize();
        state.counters["free_velocity_dofs"] = baseline.freeVelocity;
        state.counters["zero_mean_pressure_dofs"] = baseline.zeroMeanPressure;
        state.counters["velocity_degree"] = K;
        state.counters["pressure_degree"] = K - 1;
        state.counters["geometry_degree"] = curved ? 2 : 1;
        state.counters["quadrature_order"] = order;
        state.counters["quadrature_points_per_cell"] =
          QF::PolytopeQuadratureFormula::get(order, geometry).getSize();
        state.counters["velocity_nnz"] = a.getOperator().nonZeros();
        state.counters["divergence_nnz"] = b.getOperator().nonZeros();
        state.counters["pressure_mass_nnz"] = mass.getOperator().nonZeros();
      }

      /// @brief Registers three hierarchy levels on six geometries and two maps.
      static void registerCases()
      {
        for (const auto geometry :
          {Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
            Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
            Polytope::Type::Hexahedron, Polytope::Type::Wedge})
        {
          for (const bool curved : {false, true})
          {
            const std::string name = std::string("Stokes/PressureOracle/") +
              (curved ? "CurvedQ2/" : "Affine/") + "P" + std::to_string(K) + "P" +
              std::to_string(K - 1) + "/" +
              std::string(UniformGrid::getGeometryName(geometry));
            auto* entry = benchmark::RegisterBenchmark(
              name.c_str(), [geometry, curved](benchmark::State& state) {
                run(state, geometry, curved);
              });
            const std::array<int, 3> levels =
              K == 2 ? std::array<int, 3>{2, 3, 5} : std::array<int, 3>{3, 4, 5};
            for (const int n : levels)
            {
              for (const int order : {16, 20})
                entry->Args({n, order});
            }
            entry->UseRealTime();
          }
        }
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  ::testing::InitGoogleTest(&argc, argv);
  PressureOracle<2>::registerCases();
  PressureOracle<3>::registerCases();
  PressureOracle<4>::registerCases();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("scope", "Pressure-spectrum oracle; FE assembly excluded");
  benchmark::AddCustomContext("compiler", __VERSION__);
  const auto count = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return !count || pressureOracleFailed || ::testing::Test::HasFailure() ? 1 : 0;
}
