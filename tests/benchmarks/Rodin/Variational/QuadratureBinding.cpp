/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Local binding of the operator terms consumed by coupled assembly.
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include "Workloads/CoupledReactionDiffusionWorkload.h"

namespace Rodin::Tests::Benchmarks
{
  bool quadratureBindingFailed = false;

  /**
   * @brief Profiles complete local preparation through public integral binding.
   *
   * ## Architecture
   *
   * Each scope binds one of the four distinct operator types in the coupled
   * reaction-diffusion workload. It includes formula selection, owned mapped
   * quadrature, basis preparation and local matrix contraction; it excludes
   * global insertion, constraints and solving. The mesh, element degrees and
   * order match the Assembly-owned consumer. Reference-table first use and
   * numerical checks are outside timing. Mesh traversal remains included.
   *
   * Every cell matrix is qualified against the existing original entry-loop
   * oracle, which shares mapped geometry but not production basis reuse. A
   * perturbed entry must fail that comparison. A separate untimed sweep after
   * measurement must reproduce every initially retained entry exactly.
   *
   * @tparam K Scalar H1 degree.
   * @tparam Scalar Real or complex field.
   * @tparam Term Mass, weighted mass, diffusion or weighted diffusion index.
   */
  template <size_t K, class Scalar, size_t Term>
  class QuadratureBinding
  {
    public:
      /// @brief Matches the complete assembly oracle's relative matrix budget.
      static constexpr Real ReferenceTolerance = 1e-11;

      /**
       * @brief Measures repeated cell binding after independent entry qualification.
       * @param state Grid points per axis and repetition protocol.
       * @param geometry Reference-cell family for the unit-box mesh.
       */
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry)
      {
        using namespace Variational;
        using Workload = CoupledReactionDiffusionWorkload<K, Scalar>;
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
        H1<K, Scalar> space(std::integral_constant<size_t, K>{}, mesh);
        TrialFunction u(space);
        TestFunction v(space);
        auto integral = [&]() {
          if constexpr (Term == 0)
            return Integral(u, v);
          else if constexpr (Term == 1)
            return Integral(Convergence::ReactionDiffusionData::Coupling * u, v);
          else if constexpr (Term == 2)
            return Integral(Grad(u), Grad(v));
          else
            return Integral(Real(2) * Grad(u), Grad(v));
        }();
        integral.setOrder(Workload::QuadratureOrder);
        QuadratureReference reference(integral);
        std::vector<Math::Matrix<Scalar>> initial;
        Real maxReferenceError = 0;
        size_t localEntries = 0;
        const auto matches = [](const auto& lhs, const auto& rhs) {
          return lhs.rows() == rhs.rows() && lhs.cols() == rhs.cols() &&
            (lhs - rhs).norm() <= ReferenceTolerance * std::max(Real(1), rhs.norm());
        };
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          integral.setPolytope(*cell);
          reference.assemble(*cell, Workload::QuadratureOrder);
          const auto& expected = reference.getOperator();
          Math::Matrix<Scalar> actual(expected.rows(), expected.cols());
          for (Eigen::Index te = 0; te < actual.rows(); ++te)
          {
            for (Eigen::Index tr = 0; tr < actual.cols(); ++tr)
              actual(te, tr) = integral.integrate(tr, te);
          }
          auto wrong = expected;
          wrong(0, 0) += std::max(Real(1), expected.norm());
          if (!matches(actual, expected) || matches(wrong, expected))
          {
            quadratureBindingFailed = true;
            state.SkipWithError("Local matrix or perturbation-rejection guard failed");
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
              const Scalar value = integral.integrate(tr, te);
              if (!std::isfinite(std::abs(value)) || value != expected(te, tr))
              {
                quadratureBindingFailed = true;
                state.SkipWithError("Repeated local binding changed an initial entry");
                return;
              }
            }
          }
        }
        const size_t d = mesh.getDimension();
        state.counters["cells"] = mesh.getPolytopeCount(d);
        state.counters["dimension"] = d;
        state.counters["degree"] = K;
        state.counters["global_dofs"] = space.getSize();
        state.counters["local_basis"] = space.getFiniteElement(d, 0).getCount();
        state.counters["local_entries"] = localEntries;
        state.counters["quadrature_order"] = Workload::QuadratureOrder;
        state.counters["quadrature_points_per_cell"] =
          QF::PolytopeQuadratureFormula::get(Workload::QuadratureOrder, geometry)
            .getSize();
        state.counters["relative_reference_error"] = maxReferenceError;
        state.counters["specialized"] = decltype(integral)::Specialized;
      }

      /// @brief Registers matched consumer geometries and three mesh sizes.
      static void registerCases()
      {
        constexpr const char* scope = Term == 0 ? "Binding/Mass"
          : Term == 1                           ? "Binding/WeightedMass"
          : Term == 2                           ? "Binding/Diffusion"
                                                : "Binding/WeightedDiffusion";
        CoupledReactionDiffusionWorkload<K, Scalar>::registerCases(
          scope, [](auto& state, auto geometry) { run(state, geometry); });
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  const auto registerDegree = []<size_t K>() {
    const auto registerScalar = []<class Scalar>() {
      QuadratureBinding<K, Scalar, 0>::registerCases();
      QuadratureBinding<K, Scalar, 1>::registerCases();
      QuadratureBinding<K, Scalar, 2>::registerCases();
      QuadratureBinding<K, Scalar, 3>::registerCases();
    };
    registerScalar.template operator()<Real>();
    registerScalar.template operator()<Complex>();
  };
  registerDegree.template operator()<1>();
  registerDegree.template operator()<2>();
  registerDegree.template operator()<3>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext("scope", "Local binding sweep; no scatter or solve");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || quadratureBindingFailed ? 1 : 0;
}
