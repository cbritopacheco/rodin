/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include "Rodin/Assembly.h"
#include "QuadratureReference.h"
#include "convergence/Stokes.h"

namespace Rodin::Tests::Benchmarks
{
  bool stokesBindingFailed = false;

  /**
   * @brief Isolates the three local operator kernels consumed by mixed Stokes assembly.
   *
   * ## Architecture
   *
   * The terms are @f$(\nabla u,\nabla v)@f$, @f$(\nabla\cdot u,q)@f$
   * and @f$(p,\nabla\cdot v)@f$, with velocity degree K and pressure degree
   * K-1. Public integral binding includes owned mapped quadrature, basis
   * transforms and eager local contraction. Mesh traversal and one entry read
   * remain timed; global insertion, constraints and solving do not.
   * Order 16 remains the default and matches the Assembly-owned Stokes workload.
   * Separately named boundary studies hold each mesh and mixed pair fixed while
   * varying only the integration order; they do not change convergence tests.
   * The original generic
   * entry loop qualifies every cell matrix before timing; an entry perturbation
   * must be rejected. A final sweep must reproduce each retained entry exactly.
   * Both mixed orientations are measured rather than inferred from transposition.
   *
   * @tparam K Velocity degree, at least two.
   * @tparam Term Zero for vector diffusion, one for divergence, two for pressure,
   * or three for explicitly unit-weighted vector diffusion.
   */
  template <size_t K, size_t Term>
  class StokesBinding
  {
    public:
      /// @brief Shared polynomial quadrature order of the complete Stokes consumer.
      static constexpr size_t QuadratureOrder = 16;

      /// @brief Relative Frobenius budget used by the original-loop assembly oracle.
      static constexpr Real ReferenceTolerance = 1e-11;

      /**
       * @brief Measures local preparation after qualifying every matrix entry.
       * @param state Grid points per axis and repetition protocol.
       * @param geometry Unit-box reference-cell family.
       * @param quadratureOrder Explicit integration degree; default sixteen.
       */
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry,
        size_t quadratureOrder = QuadratureOrder)
      {
        using namespace Variational;
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
        const size_t dim = mesh.getDimension();
        H1 velocity(std::integral_constant<size_t, K>{}, mesh, dim);
        H1 pressure(std::integral_constant<size_t, K - 1>{}, mesh);
        TrialFunction u(velocity);
        TrialFunction p(pressure);
        TestFunction v(velocity);
        TestFunction q(pressure);
        const auto integrand = [&]() {
          if constexpr (Term == 0)
            return Dot(Jacobian(u), Jacobian(v));
          else if constexpr (Term == 1)
            return Dot(Div(u), q);
          else if constexpr (Term == 2)
            return Dot(p, Div(v));
          else
            return Dot(Real(1) * Jacobian(u), Jacobian(v));
        }();
        Integral integral(integrand);
        integral.setOrder(quadratureOrder);
        QuadratureReference<decltype(integral)> reference(integrand);
        std::vector<Math::Matrix<Real>> initial;
        Real maxReferenceError = 0;
        size_t localEntries = 0;
        const auto matches = [](const auto& lhs, const auto& rhs) {
          return lhs.rows() == rhs.rows() && lhs.cols() == rhs.cols() &&
            (lhs - rhs).norm() <= ReferenceTolerance * std::max(Real(1), rhs.norm());
        };
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          integral.setPolytope(*cell);
          reference.assemble(*cell, quadratureOrder);
          const auto& expected = reference.getOperator();
          Math::Matrix<Real> actual(expected.rows(), expected.cols());
          for (Eigen::Index te = 0; te < actual.rows(); ++te)
          {
            for (Eigen::Index tr = 0; tr < actual.cols(); ++tr)
              actual(te, tr) = integral.integrate(tr, te);
          }
          auto wrong = expected;
          wrong(0, 0) += std::max(Real(1), expected.norm());
          if (!matches(actual, expected) || matches(wrong, expected))
          {
            stokesBindingFailed = true;
            state.SkipWithError("Stokes local matrix or perturbation guard failed");
            return;
          }
          maxReferenceError = std::max(maxReferenceError,
            (actual - expected).norm() / std::max(Real(1), expected.norm()));
          localEntries += actual.size();
          initial.emplace_back(std::move(actual));
        }
        for (auto _ : state)
        {
          for (auto cell = mesh.getCell(); cell; ++cell)
          {
            integral.setPolytope(*cell);
            benchmark::DoNotOptimize(integral.integrate(0, 0));
          }
          benchmark::ClobberMemory();
        }
        size_t index = 0;
        for (auto cell = mesh.getCell(); cell; ++cell, ++index)
        {
          integral.setPolytope(*cell);
          const auto& expected = initial[index];
          for (Eigen::Index te = 0; te < expected.rows(); ++te)
          {
            for (Eigen::Index tr = 0; tr < expected.cols(); ++tr)
            {
              const Real value = integral.integrate(tr, te);
              if (!std::isfinite(value) || value != expected(te, tr))
              {
                stokesBindingFailed = true;
                state.SkipWithError("Repeated Stokes binding changed an initial entry");
                return;
              }
            }
          }
        }
        state.counters["cells"] = mesh.getPolytopeCount(dim);
        state.counters["dimension"] = dim;
        state.counters["velocity_degree"] = K;
        state.counters["pressure_degree"] = K - 1;
        state.counters["velocity_dofs"] = velocity.getSize();
        state.counters["pressure_dofs"] = pressure.getSize();
        state.counters["local_rows"] = initial.front().rows();
        state.counters["local_columns"] = initial.front().cols();
        state.counters["local_entries"] = localEntries;
        state.counters["quadrature_order"] = quadratureOrder;
        const size_t samples =
          QF::PolytopeQuadratureFormula::get(quadratureOrder, geometry).getSize();
        state.counters["quadrature_points_per_cell"] = samples;
        state.counters["cell_sample_pairs"] = mesh.getPolytopeCount(dim) * samples;
        state.counters["local_entry_sample_pairs"] = localEntries * samples;
        state.counters["relative_reference_error"] = maxReferenceError;
        state.counters["specialized"] = decltype(integral)::Specialized;
      }

      /**
       * @brief Registers fixed-order mesh scaling and separate boundary studies.
       *
       * The original three mesh sizes retain degree sixteen. Boundary studies
       * use two grid points per axis and degrees fifteen through seventeen;
       * tetrahedra additionally retain degree fourteen and pyramids ten and
       * eleven. These small fixed meshes isolate local-kernel scaling rather
       * than estimate global insertion or solve costs. Entry/sample pairs are
       * a logical dense-operator work proxy, not executed arithmetic counts.
       * ScalarBlockScaling additionally varies orders zero, four and eight
       * for the two vector-diffusion kernels on a fixed two-point mesh. These
       * are discrete local-operator cost studies, not integration-exactness or
       * Stokes convergence assertions; complete consumers retain order sixteen.
       */
      static void registerCases()
      {
        using Geometry::Polytope;
        for (const auto geometry :
          {Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
            Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
            Polytope::Type::Hexahedron, Polytope::Type::Wedge})
        {
          const std::string name = std::string("Stokes/Binding/") +
            (Term == 0      ? "VectorDiffusion/"
                : Term == 1 ? "Divergence/"
                : Term == 2 ? "Pressure/"
                            : "WeightedVectorDiffusion/") +
            "H1/P" + std::to_string(K) + "P" + std::to_string(K - 1) + "/Real/" +
            std::string(Convergence::UniformGrid::getGeometryName(geometry));
          benchmark::RegisterBenchmark(
            name.c_str(), [geometry](auto& state) { run(state, geometry); })
            ->Arg(3)
            ->Arg(5)
            ->Arg(9)
            ->UseRealTime();
          const std::string boundaryName = name + "/QuadratureBoundary";
          auto* boundary = benchmark::RegisterBenchmark(boundaryName.c_str(),
            [geometry](auto& state) { run(state, geometry, state.range(1)); });
          for (const int order : {15, 16, 17})
            boundary->Args({2, order});
          if (geometry == Polytope::Type::Tetrahedron)
            boundary->Args({2, 14});
          else if (geometry == Polytope::Type::Pyramid)
          {
            for (const int order : {10, 11})
              boundary->Args({2, order});
          }
          boundary->UseRealTime();
          if constexpr (Term == 0 || Term == 3)
          {
            const std::string scalingName = name + "/ScalarBlockScaling";
            auto* scaling = benchmark::RegisterBenchmark(scalingName.c_str(),
              [geometry](auto& state) { run(state, geometry, state.range(1)); });
            for (const int order : {0, 4, 8})
              scaling->Args({2, order});
            scaling->UseRealTime();
          }
        }
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  const auto registerDegree = []<size_t K>() {
    StokesBinding<K, 0>::registerCases();
    StokesBinding<K, 1>::registerCases();
    StokesBinding<K, 2>::registerCases();
    StokesBinding<K, 3>::registerCases();
  };
  registerDegree.template operator()<2>();
  registerDegree.template operator()<3>();
  registerDegree.template operator()<4>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext("scope", "Local Stokes binding; no scatter or solve");
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || stokesBindingFailed ? 1 : 0;
}
