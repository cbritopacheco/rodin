/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Separates pyramid index, Bernstein and modal source-kernel costs.
 */
#include <benchmark/benchmark.h>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include "Rodin/Variational/H1/H1Element.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

namespace Rodin::Tests::Benchmarks
{
  bool pyramidFailed = false;

  enum class PyramidOperation
  {
    Index,
    Binomial,
    Bernstein,
    BernsteinDerivative,
    Modal,
    ModalDerivative
  };

  /**
   * @brief Attributes intrinsic real pyramid kernels without nodal contraction.
   *
   * Metadata, Bernstein factors and modal evaluation are measured separately.
   * No Vandermonde inverse, reference-table lookup, map, mesh, assembly or solve
   * is timed. Degrees one through six cover the supported performance matrix.
   * Factor and modal sweeps use orders two, four and eight; metadata sweeps use
   * their degree-dependent index sets instead of meaningless quadrature sizes.
   * Runtime argument arrays are exposed through compiler barriers before each
   * sweep so that constant loop bounds cannot turn binomial or index decoding
   * into precomputed outputs. The barrier and output-write costs remain in the
   * reported sweep time; these are not isolated instruction latencies.
   *
   * ## Independent numerical oracle
   *
   * Binomial coefficients are constructed by Pascal addition, independently of
   * the source's ordered product/division. Factors use
   * @f$B_i^n(t)=\binom ni t^i(1-t)^{n-i}@f$ and direct polynomial differentiation,
   * evaluated with powers rather than the production multiplication loops.
   * With @f$a=x/(1-z)@f$, @f$b=y/(1-z)@f$, the modal oracle is
   * @f$B_i^{K-k}(a)B_j^{K-k}(b)B_k^K(z)@f$; its gradient follows by the chain
   * rule. Mode triples are enumerated independently in layer order, not decoded
   * through the measured implementation. An absolute/relative roundoff budget
   * of @f$10^{-12}\max(1,|q_{\rm ref}|)@f$ applies to dimensionless values.
   * Metadata comparisons are exact, and a perturbed output must be rejected.
   * Real kernels also underlie complex finite-element coefficients; complex
   * expansion and nodal contraction belong to the reference-element profiler.
   *
   * @tparam K Polynomial degree.
   */
  template <size_t K>
  class Pyramid
  {
    public:
      /// @brief Roundoff allowance for independent polynomial evaluation.
      static constexpr Real OracleTolerance = 1e-12;

      /**
       * @brief Measures one source operation and checks every written output.
       * @param state Quadrature order where applicable and timing iterations.
       * @param operation Selected metadata, factor or modal scope.
       */
      static void run(benchmark::State& state, PyramidOperation operation)
      {
        using namespace Variational;
        const bool metadata =
          operation == PyramidOperation::Index || operation == PyramidOperation::Binomial;
        const bool factor = operation == PyramidOperation::Bernstein ||
          operation == PyramidOperation::BernsteinDerivative;
        const bool derivative = operation == PyramidOperation::BernsteinDerivative ||
          operation == PyramidOperation::ModalDerivative;
        const QF::QuadratureFormulaBase* formula = nullptr;
        if (!metadata)
          formula = &QF::PolytopeQuadratureFormula::get(
            state.range(0), Geometry::Polytope::Type::Pyramid);
        const size_t nqp = formula ? formula->getSize() : 0;
        std::vector<std::array<size_t, 3>> modes;
        std::vector<size_t> modeIndices;
        std::vector<std::array<size_t, 2>> factorArguments;
        for (size_t k = 0; k <= K; ++k)
        {
          for (size_t j = 0; j <= K - k; ++j)
          {
            for (size_t i = 0; i <= K - k; ++i)
              modes.push_back({i, j, k});
          }
        }
        for (size_t mode = 0; mode < modes.size(); ++mode)
          modeIndices.push_back(mode);
        for (size_t n = 0; n <= K; ++n)
        {
          for (size_t i = 0; i <= n; ++i)
            factorArguments.push_back({n, i});
        }
        std::array<std::array<Real, K + 1>, K + 1> pascal{};
        for (size_t n = 0; n <= K; ++n)
        {
          pascal[n][0] = pascal[n][n] = 1;
          for (size_t i = 1; i < n; ++i)
            pascal[n][i] = pascal[n - 1][i - 1] + pascal[n - 1][i];
        }
        const auto polynomial = [&](size_t n, size_t i, Real x) {
          return pascal[n][i] * std::pow(x, i) * std::pow(1 - x, n - i);
        };
        const auto polynomialDerivative = [&](size_t n, size_t i, Real x) {
          Real value = 0;
          if (i > 0)
            value += pascal[n][i] * Real(i) * std::pow(x, i - 1) * std::pow(1 - x, n - i);
          if (i < n)
            value -=
              pascal[n][i] * Real(n - i) * std::pow(x, i) * std::pow(1 - x, n - i - 1);
          return value;
        };
        std::vector<Real> expected;
        if (operation == PyramidOperation::Index)
        {
          for (const auto& ijk : modes)
          {
            for (size_t j = 0; j < 3; ++j)
              expected.push_back(ijk[j]);
          }
        }
        else if (operation == PyramidOperation::Binomial)
        {
          for (size_t n = 0; n <= K; ++n)
          {
            for (size_t i = 0; i <= n; ++i)
              expected.push_back(pascal[n][i]);
          }
        }
        else
        {
          for (size_t qp = 0; qp < nqp; ++qp)
          {
            const auto& r = formula->getPoint(qp);
            if (factor)
            {
              for (size_t n = 0; n <= K; ++n)
              {
                for (size_t i = 0; i <= n; ++i)
                {
                  for (size_t j = 0; j < 3; ++j)
                  {
                    expected.push_back(derivative ? polynomialDerivative(n, i, r(j))
                                                  : polynomial(n, i, r(j)));
                  }
                }
              }
            }
            else
            {
              const Real q = 1 - r.z(), a = r.x() / q, b = r.y() / q;
              for (const auto& ijk : modes)
              {
                const auto [i, j, k] = ijk;
                const size_t n = K - k;
                const Real ba = polynomial(n, i, a), bb = polynomial(n, j, b);
                const Real bz = polynomial(K, k, r.z());
                if (!derivative)
                  expected.push_back(ba * bb * bz);
                else
                {
                  const Real da = polynomialDerivative(n, i, a);
                  const Real db = polynomialDerivative(n, j, b);
                  const Real dz = polynomialDerivative(K, k, r.z());
                  expected.push_back(da * bb * bz / q);
                  expected.push_back(ba * db * bz / q);
                  expected.push_back((da * a * bb + ba * db * b) * bz / q + ba * bb * dz);
                }
              }
            }
          }
        }
        std::vector<Real> output(expected.size());
        const auto evaluate = [&]() {
          benchmark::DoNotOptimize(modeIndices.data());
          benchmark::DoNotOptimize(factorArguments.data());
          benchmark::ClobberMemory();
          size_t pos = 0;
          if (operation == PyramidOperation::Index)
          {
            for (const size_t mode : modeIndices)
            {
              size_t i, j, k;
              PyramidIndex<K>::decode(mode, i, j, k);
              output[pos++] = i;
              output[pos++] = j;
              output[pos++] = k;
            }
          }
          else if (operation == PyramidOperation::Binomial)
          {
            for (const auto& argument : factorArguments)
            {
              const auto [n, i] = argument;
              output[pos++] = BernsteinPyramid<K>::getBinomial(n, i);
            }
          }
          else
          {
            for (size_t qp = 0; qp < nqp; ++qp)
            {
              const auto& r = formula->getPoint(qp);
              if (factor)
              {
                for (const auto& argument : factorArguments)
                {
                  const auto [n, i] = argument;
                  for (size_t j = 0; j < 3; ++j)
                  {
                    output[pos++] = derivative
                      ? BernsteinPyramid<K>::getDerivative(n, i, r(j))
                      : BernsteinPyramid<K>::getBasis(n, i, r(j));
                  }
                }
              }
              else
              {
                for (const size_t mode : modeIndices)
                {
                  if (!derivative)
                    output[pos++] = PyramidModal<K>::getBasis(mode, r);
                  else
                  {
                    for (size_t j = 0; j < 3; ++j)
                      output[pos++] = PyramidModal<K>::getDerivative(mode, j, r);
                  }
                }
              }
            }
          }
        };
        const auto valid = [&]() {
          for (size_t i = 0; i < expected.size(); ++i)
          {
            if (!std::isfinite(output[i]) || (metadata && output[i] != expected[i]) ||
              std::abs(output[i] - expected[i]) >
                OracleTolerance * std::max(Real(1), std::abs(expected[i])))
              return false;
          }
          return true;
        };
        const auto fail = [&](const char* message) {
          pyramidFailed = true;
          state.SkipWithError(message);
        };
        evaluate();
        if (!valid())
        {
          fail("Pyramid kernel disagrees with independent polynomial oracle");
          return;
        }
        output.front() += 1;
        if (valid())
        {
          fail("Pyramid oracle accepted a perturbed output");
          return;
        }
        evaluate();
        for (auto _ : state)
        {
          evaluate();
          benchmark::DoNotOptimize(output.data());
          benchmark::ClobberMemory();
        }
        if (!valid())
          fail("Timed pyramid kernel changed numerical values");
        state.counters["degree"] = K;
        state.counters["modes"] = modes.size();
        state.counters["quadrature_points"] = nqp;
        state.counters["output_entries"] = output.size();
      }

      /// @brief Registers metadata and three quadrature workloads separately.
      static void registerCases()
      {
        for (const auto& [operation, label] :
          {std::pair{PyramidOperation::Index, "IndexDecode"},
            {PyramidOperation::Binomial, "Binomial"},
            {PyramidOperation::Bernstein, "BernsteinValues"},
            {PyramidOperation::BernsteinDerivative, "BernsteinDerivatives"},
            {PyramidOperation::Modal, "ModalValues"},
            {PyramidOperation::ModalDerivative, "ModalDerivatives"}})
        {
          const std::string name =
            "H1/Pyramid/" + std::string(label) + "/P" + std::to_string(K);
          auto* registration = benchmark::RegisterBenchmark(
            name.c_str(), [operation](auto& state) { run(state, operation); });
          if (operation != PyramidOperation::Index &&
            operation != PyramidOperation::Binomial)
            registration->Arg(2)->Arg(4)->Arg(8);
          registration->UseRealTime();
        }
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  Pyramid<1>::registerCases();
  Pyramid<2>::registerCases();
  Pyramid<3>::registerCases();
  Pyramid<4>::registerCases();
  Pyramid<5>::registerCases();
  Pyramid<6>::registerCases();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || pyramidFailed ? 1 : 0;
}
