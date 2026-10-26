/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Full Eigen mass/reaction operator and projection-load assembly.
 */
#include "MassAssembly.h"

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  registerMass<Rodin::Real, false>();
  registerMass<Rodin::Real, true>();
  registerMass<Rodin::Complex, false>();
  registerMass<Rodin::Complex, true>();
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
  benchmark::AddCustomContext(
    "scope", "Warmed operator or projection-load assembly; no solve");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return !cases || massCheckFailed ? 1 : 0;
}
