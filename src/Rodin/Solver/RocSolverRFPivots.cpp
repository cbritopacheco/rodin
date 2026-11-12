/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <hip/hip_runtime.h>

#include "RocSolverRFPivots.h"

namespace Rodin::Solver::Internal
{
  /**
   * @brief Combines invalid pivots into one device flag.
   * @param values Gathered device pivots.
   * @param count Number of pivots.
   * @param flag Device flag initialized to zero before this launch.
   */
  __global__ void checkRocSolverRFPivots(const double* values, size_t count, int* flag)
  {
    const size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < count && !RocSolverRFPivots::valid(values[i]))
      atomicExch(flag, 1);
  }

  hipError_t RocSolverRFPivots::check(
    const double* values, size_t count, int* flag, bool& result, hipStream_t stream)
  {
    result = false;
    hipError_t status = hipMemsetAsync(flag, 0, sizeof(int), stream);
    if (status != hipSuccess)
      return status;
    if (count != 0)
    {
      // Four CDNA waves or eight RDNA waves per block; no shared storage.
      constexpr unsigned Threads = 256;
      const unsigned blocks = static_cast<unsigned>((count - 1) / Threads + 1);
      hipLaunchKernelGGL(checkRocSolverRFPivots, dim3(blocks), dim3(Threads), 0, stream,
        values, count, flag);
      status = hipGetLastError();
      if (status != hipSuccess)
        return status;
    }
    int invalid = 0;
    status =
      hipMemcpyAsync(&invalid, flag, sizeof(invalid), hipMemcpyDeviceToHost, stream);
    // Complete submitted work even on a copy failure, before stack storage dies.
    const hipError_t completion = hipStreamSynchronize(stream);
    if (status != hipSuccess)
      return status;
    if (completion != hipSuccess)
      return completion;
    result = invalid == 0;
    return hipSuccess;
  }
}
