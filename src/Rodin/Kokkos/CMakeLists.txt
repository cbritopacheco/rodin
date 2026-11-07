#[[
         Copyright Carlos BRITO PACHECO 2026.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

add_library(RodinKokkos INTERFACE)
add_library(Rodin::Kokkos ALIAS RodinKokkos)
set_target_properties(RodinKokkos PROPERTIES EXPORT_NAME Kokkos)

target_compile_features(RodinKokkos INTERFACE cxx_std_20)

target_link_libraries(RodinKokkos
  INTERFACE
    Rodin::Geometry
    Rodin::Context
    Rodin::Variational
    Rodin::Assembly
    Kokkos::kokkos
)

# IGC gives a thread 128 registers by default on the Arc Pro B70, and the
# bilinear integration kernels spill hundreds of them, which made them four
# times slower than on an RTX 2080 Ti. The automatic mode gives 256 registers
# to the kernels that would spill and leaves the others at 128. The option goes
# to ocloc, which compiles the device code ahead of time; SHELL keeps CMake
# from merging the repeated -Xsycl-target-backend with the one of Kokkos.
if (Kokkos_ENABLE_SYCL AND "INTEL_BMG" IN_LIST Kokkos_ARCH)
  target_compile_options(RodinKokkos INTERFACE
    "$<$<COMPILE_LANGUAGE:CXX>:SHELL:-Xsycl-target-backend \"-options -ze-intel-enable-auto-large-GRF-mode\">")
endif()

target_sources(RodinKokkos INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ForwardDecls.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Context.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Context/CUDA.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Context/HIP.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Context/SYCL.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Context/Kokkos.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Geometry.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/LinearSystem.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SparseMatrix.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Assembly.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Assembly/Kokkos.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Assembly/KokkosIteration.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Utility.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Geometry/Connectivity.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Geometry/Mesh.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Geometry/MeshView.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/DOFs.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/DeviceView.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/H1.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/HyperElasticLaw.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/Integral.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/InternalVirtualWork.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/LinearElasticityIntegral.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/P0g.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/P1.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/Potential.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/Projection.h>
)

install(TARGETS RodinKokkos EXPORT RodinTargets)

if (RODIN_USE_PETSC)
  add_library(RodinKokkosPETSc INTERFACE)
  add_library(Rodin::KokkosPETSc ALIAS RodinKokkosPETSc)
  set_target_properties(RodinKokkosPETSc PROPERTIES EXPORT_NAME KokkosPETSc)

  target_link_libraries(RodinKokkosPETSc
    INTERFACE
      Rodin::Kokkos
      Rodin::PETSc
  )

  target_sources(RodinKokkosPETSc INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/PETSc.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Assembly/PETSc.h>
  )

  install(TARGETS RodinKokkosPETSc EXPORT RodinTargets)

  if (RODIN_USE_MPI)
    add_library(RodinKokkosMPI INTERFACE)
    add_library(Rodin::KokkosMPI ALIAS RodinKokkosMPI)
    set_target_properties(RodinKokkosMPI PROPERTIES EXPORT_NAME KokkosMPI)

    target_link_libraries(RodinKokkosMPI
      INTERFACE
        Rodin::KokkosPETSc
        Rodin::MPI
    )

    target_sources(RodinKokkosMPI INTERFACE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MPI.h>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MPI/Geometry/Mesh.h>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MPI/Variational/P1.h>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MPI/Variational/P0g.h>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MPI/Variational/H1.h>
    )

    install(TARGETS RodinKokkosMPI EXPORT RodinTargets)
  endif()
endif()

if (RODIN_USE_KOKKOSKERNELS)
  add_library(RodinKokkosKernels INTERFACE)
  add_library(Rodin::KokkosKernels ALIAS RodinKokkosKernels)
  set_target_properties(RodinKokkosKernels PROPERTIES EXPORT_NAME KokkosKernels)

  target_link_libraries(RodinKokkosKernels
    INTERFACE
      Rodin::Kokkos
      Rodin::Solver
      Kokkos::kokkoskernels
  )

  target_sources(RodinKokkosKernels INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solver.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solver/Krylov.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solver/CG.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solver/GMRES.h>
  )

  install(TARGETS RodinKokkosKernels EXPORT RodinTargets)
endif()

if (RODIN_USE_HYPRE)
  add_library(RodinKokkosHypre INTERFACE)
  add_library(Rodin::KokkosHypre ALIAS RodinKokkosHypre)
  set_target_properties(RodinKokkosHypre PROPERTIES EXPORT_NAME KokkosHypre)

  target_link_libraries(RodinKokkosHypre
    INTERFACE
      Rodin::Kokkos
      Rodin::Solver
      HYPRE::HYPRE
  )

  target_sources(RodinKokkosHypre INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solver/BoomerAMG.h>
  )

  install(TARGETS RodinKokkosHypre EXPORT RodinTargets)
endif()

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Kokkos
)
