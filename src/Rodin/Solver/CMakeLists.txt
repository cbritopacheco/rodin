#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinSolver_HEADERS
  Info.h
  LinearSolver.h
  NewtonSolver.h
  CG.h
  SparseLU.h
)

# The Eigen solvers used by device code are instantiated here by the host
# compiler. A translation unit compiled by nvcc defines EIGEN_GPUCC, which
# turns off Eigen's vectorization in host code too, so instantiating them
# there would run the factorizations and the Krylov iterations on scalar code.
set(RodinSolver_SRCS
  CG.cpp
  SparseLU.cpp
)

add_library(RodinSolver ${RodinSolver_SRCS} ${RodinSolver_HEADERS})
add_library(Rodin::Solver ALIAS RodinSolver)
set_target_properties(RodinSolver PROPERTIES EXPORT_NAME Solver)

# Consumers use the instantiations above instead of instantiating their own.
target_compile_definitions(RodinSolver PUBLIC RODIN_SOLVER_EXTERN_TEMPLATES)

# ---- OpenMP ----------------------------------------------------------------
if (RODIN_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(RodinSolver PUBLIC OpenMP::OpenMP_CXX)
endif()

# ---- SuiteSparse -----------------------------------------------------------
if (RODIN_USE_UMFPACK)
  target_link_libraries(RodinSolver PUBLIC SuiteSparse::UMFPACK)
endif()

if (RODIN_USE_PARU)
  target_link_libraries(RodinSolver PUBLIC SuiteSparse::ParU)
endif()

# ---- MUMPS -----------------------------------------------------------------
if (RODIN_USE_MUMPS)
  target_link_libraries(RodinSolver PUBLIC MUMPS::MUMPS)
endif()

# ---- cuDSS -----------------------------------------------------------------
if (RODIN_USE_CUDSS)
  target_link_libraries(RodinSolver PUBLIC cudss CUDA::cudart)
endif()

# ---- rocSOLVER -------------------------------------------------------------
if (RODIN_USE_ROCSOLVER)
  target_sources(RodinSolver PRIVATE RocSolverRFPivots.cpp)
  set_source_files_properties(RocSolverRFPivots.cpp PROPERTIES LANGUAGE HIP)
  set_target_properties(RodinSolver PROPERTIES HIP_STANDARD 17 HIP_STANDARD_REQUIRED YES)
  target_link_libraries(RodinSolver PUBLIC roc::rocsolver roc::rocsparse roc::rocblas hip::host)
endif()

if (RODIN_USE_SPQR)
  target_link_libraries(RodinSolver PUBLIC SuiteSparse::SPQR)
endif()

if (RODIN_USE_CHOLMOD)
  target_link_libraries(RodinSolver PUBLIC SuiteSparse::CHOLMOD)
endif()

if (RODIN_USE_APPLE_ACCELERATE)
  # Header only, exposed via interface sources
  target_sources(RodinSolver PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/AppleAccelerate.h>
  )
endif()

# ---- Link targets ----------------------------------------------------------
target_link_libraries(RodinSolver PUBLIC Rodin::Rodin Rodin::Alert Eigen3::Eigen)

target_include_directories(RodinSolver
  PUBLIC
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# ---- Install target --------------------------------------------------------
install(
  TARGETS RodinSolver
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Solver
)
