/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Controlled nodal derivative contractions with and without modal reuse.
 */
#include <benchmark/benchmark.h>
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "Rodin/Variational/H1/H1Element.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

namespace Rodin::Tests::Benchmarks
{
  bool pyramidNodalFailed = false;

  enum class NodalProtocol
  {
    Public,
    RepeatedModes,
    SharedModes
  };

  /**
   * @brief Measures a caller-side reuse experiment without changing the library.
   *
   * ## Operation
   *
   * Nodal derivatives are @f$\partial_i\phi_a(\widehat x_q)=
   * \sum_b(V^{-1})_{ba}\partial_i\psi_b(\widehat x_q)@f$.
   * Public calls the existing element derivative. RepeatedModes spells the
   * same contraction with production modal derivatives. SharedModes computes
   * each modal derivative once per point and reuses it across nodal functions.
   * Every nodal sum retains ascending mode order; no Eigen reduction or
   * reassociation is substituted. Nodal sums use real arithmetic before
   * promotion to real or complex storage, matching the public element's
   * real-valued reference basis contract.
   *
   * ## Measurement boundary
   *
   * Vandermonde initialization, quadrature lookup, output allocation and
   * scratch allocation are untimed. Modal evaluation and nodal contractions
   * are timed, including scratch writes for SharedModes. The scratch holds
   * only one point's modal derivatives and is not a new production cache.
   * Reference-table lookup, physical mapping, assembly and solving are absent.
   * Public and RepeatedModes are separate compiler/calling-protocol controls;
   * a gain in this experiment is not yet a library or assembly speedup.
   *
   * ## Numerical guards
   *
   * Every quadrature output must exactly match the public derivative before
   * and after timing. Independent DOF-functional reconstruction recovers the
   * gradient of @f$\eta(1+x+2y+3z)@f$, with @f$\eta=1@f$ for real fields
   * and @f$\eta=1+i/2@f$ for complex fields, to an absolute dimensionless
   * roundoff budget of @f$10^{-8}@f$ at interior samples. Perturbation must be rejected.
   *
   * @tparam K Field factor degree.
   * @tparam Scalar Real or complex storage type.
   * @tparam Protocol Selected derivative-evaluation protocol.
   */
  template <size_t K, class Scalar, NodalProtocol Protocol>
  class PyramidNodal
  {
    public:
      /**
       * @brief Measures one quadrature sweep of all nodal derivatives.
       * @param state Requested quadrature order and iteration state.
       */
      static void run(benchmark::State& state)
      {
        using namespace Variational;
        constexpr Real ReconstructionTolerance = 1e-8;
        constexpr size_t Count = PyramidIndex<K>::Count;
        const H1Element<K, Scalar> element(Geometry::Polytope::Type::Pyramid);
        const auto& formula = QF::PolytopeQuadratureFormula::get(
          state.range(0), Geometry::Polytope::Type::Pyramid);
        const auto& inverse = VandermondePyramid<K>::getInverse();
        std::vector<Scalar> expected(formula.getSize() * Count * 3),
          output(expected.size());
        std::array<std::array<Real, 3>, Count> modal;
        std::array<Scalar, Count> coefficients;
        const Scalar phase = [] {
          if constexpr (std::is_same_v<Scalar, Complex>)
            return Scalar(1, 0.5);
          else
            return Scalar(1);
        }();
        for (size_t a = 0; a < Count; ++a)
        {
          coefficients[a] = element.getLinearForm(a)([phase](const auto& point) {
            return phase * Scalar(1 + point(0) + 2 * point(1) + 3 * point(2));
          });
        }
        for (size_t q = 0; q < formula.getSize(); ++q)
        {
          for (size_t a = 0; a < Count; ++a)
          {
            const auto& basis = element.getBasis(a);
            for (size_t i = 0; i < 3; ++i)
            {
              expected[(q * Count + a) * 3 + i] =
                basis.template getDerivative<1>(i)(formula.getPoint(q));
            }
          }
        }
        const auto evaluate = [&] {
          for (size_t q = 0; q < formula.getSize(); ++q)
          {
            const auto& point = formula.getPoint(q);
            if constexpr (Protocol == NodalProtocol::SharedModes)
            {
              for (size_t b = 0; b < Count; ++b)
              {
                for (size_t i = 0; i < 3; ++i)
                  modal[b][i] = PyramidModal<K>::getDerivative(b, i, point);
              }
            }
            for (size_t a = 0; a < Count; ++a)
            {
              for (size_t i = 0; i < 3; ++i)
              {
                Scalar value = Scalar(0);
                if constexpr (Protocol == NodalProtocol::Public)
                  value = element.getBasis(a).template getDerivative<1>(i)(point);
                else
                {
                  Real sum = 0;
                  for (size_t b = 0; b < Count; ++b)
                  {
                    const Real derivative = [&] {
                      if constexpr (Protocol == NodalProtocol::SharedModes)
                        return modal[b][i];
                      else
                        return PyramidModal<K>::getDerivative(b, i, point);
                    }();
                    sum += inverse(b, a) * derivative;
                  }
                  value = Scalar(sum);
                }
                output[(q * Count + a) * 3 + i] = value;
              }
            }
          }
        };
        const auto valid = [&] {
          for (size_t j = 0; j < output.size(); ++j)
          {
            if (!std::isfinite(std::abs(output[j])) || output[j] != expected[j])
              return false;
          }
          for (size_t q = 0; q < formula.getSize(); ++q)
          {
            for (size_t i = 0; i < 3; ++i)
            {
              Scalar derivative = Scalar(0);
              for (size_t a = 0; a < Count; ++a)
                derivative += coefficients[a] * output[(q * Count + a) * 3 + i];
              if (std::abs(derivative - phase * Scalar(i + 1)) > ReconstructionTolerance)
                return false;
            }
          }
          return true;
        };
        evaluate();
        if (!valid())
        {
          Real entryError = 0, reconstructionError = 0;
          for (size_t j = 0; j < output.size(); ++j)
            entryError = std::max(entryError, std::abs(output[j] - expected[j]));
          for (size_t q = 0; q < formula.getSize(); ++q)
          {
            for (size_t i = 0; i < 3; ++i)
            {
              Scalar derivative = Scalar(0);
              for (size_t a = 0; a < Count; ++a)
                derivative += coefficients[a] * output[(q * Count + a) * 3 + i];
              reconstructionError = std::max(
                reconstructionError, std::abs(derivative - phase * Scalar(i + 1)));
            }
          }
          std::ostringstream diagnostic;
          diagnostic << std::scientific << std::setprecision(17)
                     << "Nodal derivative oracle failed: entry error=" << entryError
                     << ", reconstruction error=" << reconstructionError;
          const auto message = diagnostic.str();
          pyramidNodalFailed = true;
          state.SkipWithError(message.c_str());
          return;
        }
        output.front() += Scalar(1);
        if (valid())
        {
          pyramidNodalFailed = true;
          state.SkipWithError("Nodal derivative oracle accepted a perturbation");
          return;
        }
        evaluate();
        for (auto _ : state)
        {
          benchmark::DoNotOptimize(&inverse);
          benchmark::ClobberMemory();
          evaluate();
          benchmark::DoNotOptimize(output.data());
          benchmark::ClobberMemory();
        }
        if (!valid())
        {
          pyramidNodalFailed = true;
          state.SkipWithError("Repeated nodal derivatives changed qualified values");
        }
        state.counters["factor_degree"] = K;
        state.counters["basis_functions"] = Count;
        state.counters["quadrature_points"] = formula.getSize();
        state.counters["output_entries"] = output.size();
        state.counters["modal_derivative_evaluations"] = formula.getSize() * 3 * Count *
          (Protocol == NodalProtocol::SharedModes ? 1 : Count);
        state.counters["nodal_contraction_products"] =
          formula.getSize() * 3 * Count * Count;
        state.counters["scratch_bytes"] =
          Protocol == NodalProtocol::SharedModes ? sizeof(modal) : 0;
      }

      /// @brief Registers four orders for this scalar and evaluation protocol.
      static void registerCases()
      {
        const char* protocol = Protocol == NodalProtocol::Public ? "Public"
          : Protocol == NodalProtocol::RepeatedModes             ? "RepeatedModes"
                                                                 : "SharedModes";
        const std::string name = std::string("H1/Pyramid/NodalDerivatives/") + protocol +
          "/P" + std::to_string(K) +
          (std::is_same_v<Scalar, Complex> ? "/Complex" : "/Real");
        benchmark::RegisterBenchmark(name.c_str(), &run)->Arg(2)->Arg(4)->Arg(8)->Arg(12);
      }
  };

  template <size_t K, class Scalar>
  void registerPyramidNodal()
  {
    PyramidNodal<K, Scalar, NodalProtocol::Public>::registerCases();
    PyramidNodal<K, Scalar, NodalProtocol::RepeatedModes>::registerCases();
    PyramidNodal<K, Scalar, NodalProtocol::SharedModes>::registerCases();
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  registerPyramidNodal<1, Real>();
  registerPyramidNodal<1, Complex>();
  registerPyramidNodal<2, Real>();
  registerPyramidNodal<2, Complex>();
  registerPyramidNodal<3, Real>();
  registerPyramidNodal<3, Complex>();
  registerPyramidNodal<4, Real>();
  registerPyramidNodal<4, Complex>();
  registerPyramidNodal<5, Real>();
  registerPyramidNodal<5, Complex>();
  registerPyramidNodal<6, Real>();
  registerPyramidNodal<6, Complex>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext(
    "scope", "Experimental nodal derivative sweep; no source repair");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || pyramidNodalFailed ? 1 : 0;
}
