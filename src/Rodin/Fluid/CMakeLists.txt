#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinFluid_HEADERS
  ForwardDecls.h
  Local/FlowPoint.h
  Local/Input.h
  Local/Tags.h
  Fields/StrainRate.h
  Fields/ShearRate.h
  Fields/Vorticity.h
  Fields/CauchyStress.h
  Constitutive/RheologyLaw.h
  Constitutive/Newtonian.h
  Constitutive/PowerLaw.h
  Constitutive/CarreauYasuda.h
  Constitutive/Bingham.h
  Constitutive/HerschelBulkley.h
  Integrators/ViscousForce.h
  Integrators/ViscousTangent.h
)

# Header-only / interface library
add_library(RodinFluid INTERFACE)
add_library(Rodin::Fluid ALIAS RodinFluid)
set_target_properties(RodinFluid PROPERTIES EXPORT_NAME Fluid)

# ---- Includes / defs -------------------------------------------------------
target_include_directories(RodinFluid
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinFluid
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# ---- Link targets ----------------------------------------------------------
target_link_libraries(RodinFluid
  INTERFACE
    Rodin::Geometry
    Rodin::Variational
)

# Expose headers for IDEs in build tree
foreach(_h IN LISTS RodinFluid_HEADERS)
  target_sources(RodinFluid INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

# ---- Install target --------------------------------------------------------
install(
  TARGETS RodinFluid
  EXPORT RodinTargets
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Fluid
)
