/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file TidyTU.cpp
 * @brief Translation unit aggregating the public headers for clang-tidy.
 *
 * The library is header-heavy, so per-file analysis of .cpp sources would
 * miss most of the code. This TU includes the public module headers so a
 * single clang-tidy invocation (with HeaderFilterRegex from .clang-tidy)
 * checks the whole header surface. Used by dev/style_lint and the Style CI
 * workflow; not part of the library build.
 *
 * Optional modules are included under the same feature definitions as the
 * library. The RodinTidyTU CMake target supplies their dependency include
 * paths whenever they are enabled.
 */
#include <Rodin/Configure.h>

#include <Rodin/Types.h>
#include <Rodin/Alert.h>
#include <Rodin/Math.h>
#include <Rodin/Geometry.h>
#include <Rodin/QF.h>
#include <Rodin/FormLanguage.h>
#include <Rodin/Variational.h>
#include <Rodin/Assembly.h>
#include <Rodin/Solver.h>
#include <Rodin/Solid.h>
#include <Rodin/Heart.h>
#include <Rodin/IO.h>
#include <Rodin/Utility.h>
#include <Rodin/Serialization/Export.h>
#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/Eikonal/FMM.h>
#include <Rodin/Adaptation/WNGIR.h>
#include <Rodin/MMG.h>

#ifdef RODIN_USE_MPI
#include <Rodin/MPI.h>
#endif

#ifdef RODIN_USE_PETSC
#include <Rodin/PETSc.h>
#endif

#ifdef RODIN_USE_SCOTCH
#include <Rodin/Scotch/MeshPartitioner.h>
#endif

#ifdef RODIN_USE_APPLE_ACCELERATE
#include <Rodin/Solver/AppleAccelerate.h>
#endif

int main()
{
  return 0;
}
