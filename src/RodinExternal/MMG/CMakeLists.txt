#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

include (ExternalProject)

# ---- Add main targets ------------------------------------------------------
set(RodinExternalMMG_HEADERS
  Cast.h
  Utility.h
  Configure.h.in
  ISCDProcess.h
  ForwardDecls.h
  Mesh.h
  Solution.h
  ScalarSolution.h
  VectorSolution.h
# ---- mmg2d
  MMG2D.h
  Mesh2D.h
  Advect2D.h
  Distancer2D.h
  MeshAdaptor2D.h
  MeshOptimizer2D.h
  ScalarSolution2D.h
  VectorSolution2D.h
  ImplicitDomainMesher2D.h
# ---- mmg3d
  Mesh3D.h
  ScalarSolution3D.h
# ---- mmgs
  MMGS.h
  MeshS.h
  ScalarSolutionS.h)

set(RodinExternalMMG_SRCS
  Cast.cpp
  Utility.cpp
  ISCDProcess.cpp
  ScalarSolution.cpp
# ---- mmg2d
  Mesh2D.cpp
  Advect2D.cpp
  Distancer2D.cpp
  MeshAdaptor2D.cpp
  MeshOptimizer2D.cpp
  ScalarSolution2D.cpp
  VectorSolution2D.cpp
  ImplicitDomainMesher2D.cpp
# ---- mmg3d
  Mesh3D.cpp
  ScalarSolution3D.cpp
# ---- mmgs
  MeshS.cpp
  AdvectS.cpp
  DistancerS.cpp
  ScalarSolutionS.cpp
  ImplicitDomainMesherS.cpp
  )

if (${CMAKE_BUILD_TYPE} MATCHES Release)
  rodin_add_compile_options(LANG C CXX OPTIONS -w -O3)
elseif (${CMAKE_BUILD_TYPE} MATCHES Debug)
  rodin_add_compile_options(
    LANG C CXX
    OPTIONS
    -g
    -O0
    -Wall
    -Wextra
    -Wpedantic
    -Wno-sign-compare
    -Wno-old-style-cast
    -Wno-unused-variable
    -Wno-unused-function
    -Wno-unused-parameter
    -Wno-double-promotion
    -Wno-missing-prototypes
    -Wno-missing-declarations
    -Wno-zero-as-null-pointer-constant
    -Wno-error=extra
    -Wno-error=pedantic)
endif()

add_library(RodinExternalMMG
  ${RodinExternalMMG_HEADERS} ${RodinExternalMMG_SRCS})
add_library(Rodin::External::MMG ALIAS RodinExternalMMG)

# ---- MMG ----
message(STATUS "Configuring mmg ...")
add_subdirectory(
  ${PROJECT_SOURCE_DIR}/third-party/mmg
  ${PROJECT_BINARY_DIR}/third-party/mmg)
# The previous call will readily add the mmg targets that we need to link the
# MMG wrapper. Check MmgTargets.cmake for a list of all the targets available.
# At the time of writing we have:
# - Mmg::libmmg_a
# - Mmg::libmmg2d_a
# - Mmg::libmmg3d_a
# - Mmg::libmmgs_a

# ---- ISCD::Mshdist ----
ExternalProject_Add(RodinExternalISCDMshdist
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/third-party/ISCD/Mshdist
  BINARY_DIR ${PROJECT_BINARY_DIR}/third-party/ISCD/Mshdist
  CMAKE_ARGS -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  INSTALL_COMMAND ""
  COMMENT "ISCD::Mshdist")
add_dependencies(RodinExternalMMG RodinExternalISCDMshdist)

# ---- ISCD::Advection ----
ExternalProject_Add(RodinExternalISCDAdvection
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/third-party/ISCD/Advection
  BINARY_DIR ${PROJECT_BINARY_DIR}/third-party/ISCD/Advection
  CMAKE_ARGS -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  INSTALL_COMMAND ""
  COMMENT "ISCD::Advection")
add_dependencies(RodinExternalMMG RodinExternalISCDAdvection)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(RODIN_MMG_VERBOSITY_LEVEL 1)
else()
  set(RODIN_MMG_VERBOSITY_LEVEL -1)
endif()

set(RODIN_ISCD_MSHDIST_EXECUTABLE
  ${PROJECT_BINARY_DIR}/third-party/ISCD/Mshdist/mshdist)
set(RODIN_ISCD_ADVECTION_EXECUTABLE
  ${PROJECT_BINARY_DIR}/third-party/ISCD/Advection/Advection)
configure_file(Configure.h.in Configure.h @ONLY)

target_include_directories(RodinExternalMMG
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  ${Boost_INCLUDE_DIRS}
  PUBLIC
  ${PROJECT_SOURCE_DIR}/third-party/mmg/src
  ${PROJECT_SOURCE_DIR}/third-party/mmg/include
  )
target_link_libraries(RodinExternalMMG
  PUBLIC
  Mmg::libmmg_a
  Mmg::libmmg2d_a
  Mmg::libmmg3d_a
  Mmg::libmmgs_a
  RodinAlert
  RodinMesh
  RodinVariational
  ${Boost_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT})

