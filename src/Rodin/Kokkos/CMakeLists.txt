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

target_sources(RodinKokkos INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ForwardDecls.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Context.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Context/CUDA.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Context/Kokkos.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Geometry.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Assembly.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Assembly/Kokkos.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Assembly/KokkosIteration.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Utility.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Geometry/Connectivity.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Geometry/Mesh.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Geometry/MeshView.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/DeviceView.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/Integral.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Variational/P1.h>
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
endif()

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Kokkos
)
