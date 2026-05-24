#[[
         Copyright Carlos BRITO PACHECO 2021 - 2026.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinAdaptation_HEADERS
  TargetMatrixOptimization.h
  TargetMatrixOptimization/Metrics.h
  TargetMatrixOptimization/Terms.h
  TargetMatrixOptimization/IsoparametricGeometry.h
  TargetMatrixOptimization/IsoparametricTMOPSolver.h
)

add_library(RodinAdaptation INTERFACE)
add_library(Rodin::Adaptation ALIAS RodinAdaptation)
set_target_properties(RodinAdaptation PROPERTIES EXPORT_NAME Adaptation)

target_include_directories(RodinAdaptation
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinAdaptation
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

target_link_libraries(RodinAdaptation
  INTERFACE
    Rodin::Geometry
    Rodin::Math
    Rodin::Variational
)

foreach(_h IN LISTS RodinAdaptation_HEADERS)
  target_sources(RodinAdaptation INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

install(
  TARGETS RodinAdaptation
  EXPORT RodinTargets
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Adaptation
)
