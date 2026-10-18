/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "MassBinding.h"

/**
 * @brief Registers both scalar-weighted expression orders without rewriting them.
 *
 * The shared fixture independently qualifies each operator. The dispatch
 * counter is an implementation observation, not a required truth value for
 * right multiplication. These optional measurements are not CTest gates.
 */
int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  const auto registerDegree = []<size_t K>() {
    const auto registerScalar = []<class Scalar>() {
      constexpr const char* Prefix = "H1/MassOperandOrder/";
      H1MassBinding<K, Scalar, MassRange::Scalar, true>::registerCases(Prefix);
      H1MassBinding<K, Scalar, MassRange::Scalar, true, true>::registerCases(Prefix);
      H1MassBinding<K, Scalar, MassRange::Vector, true>::registerCases(Prefix);
      H1MassBinding<K, Scalar, MassRange::Vector, true, true>::registerCases(Prefix);
      H1MassBinding<K, Scalar, MassRange::Matrix, true>::registerCases(Prefix);
      H1MassBinding<K, Scalar, MassRange::Matrix, true, true>::registerCases(Prefix);
    };
    registerScalar.template operator()<Real>();
    registerScalar.template operator()<Complex>();
  };
  registerDegree.template operator()<1>();
  registerDegree.template operator()<2>();
  registerDegree.template operator()<3>();
  registerDegree.template operator()<4>();
  registerDegree.template operator()<5>();
  registerDegree.template operator()<6>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext(
    "scope", "One-cell H1 scalar-weight operand orders; no scatter or solve");
  benchmark::AddCustomContext(
    "dispatch", "Observed specialization, not a mathematical work counter");
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || massBindingFailed ? 1 : 0;
}
