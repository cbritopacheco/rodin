/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "MassBinding.h"

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  const auto registerDegree = []<size_t K>() {
    const auto registerScalar = []<class Scalar>() {
      H1MassBinding<K, Scalar, MassRange::Scalar, false>::registerCases();
      H1MassBinding<K, Scalar, MassRange::Scalar, true>::registerCases();
      H1MassBinding<K, Scalar, MassRange::Vector, false>::registerCases();
      H1MassBinding<K, Scalar, MassRange::Vector, true>::registerCases();
      H1MassBinding<K, Scalar, MassRange::Matrix, false>::registerCases();
      H1MassBinding<K, Scalar, MassRange::Matrix, true>::registerCases();
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
    "scope", "One-cell H1 mass binding; no traversal, scatter or solve");
  benchmark::AddCustomContext(
    "work_proxies", "Logical entry/sample counts, not executed FLOPs");
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || massBindingFailed ? 1 : 0;
}
