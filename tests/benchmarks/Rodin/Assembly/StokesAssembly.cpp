/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Workloads/StokesWorkload.h"
#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  bool failed = false;
  StokesWorkload<2>::registerCases<false>(failed);
  StokesWorkload<3>::registerCases<false>(failed);
  StokesWorkload<4>::registerCases<false>(failed);
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
  benchmark::AddCustomContext(
    "timed_operation", "Public mixed Stokes assembly, constraints and completion");
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || failed ? 1 : 0;
}
