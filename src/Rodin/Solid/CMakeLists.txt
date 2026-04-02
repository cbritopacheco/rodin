#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinSolid_HEADERS
  ForwardDecls.h
  Kinematics/KinematicState.h
  Kinematics/Invariants.h
  Constitutive/HyperElasticLaw.h
  Constitutive/NeoHookean.h
  Constitutive/SaintVenantKirchhoff.h
  Constitutive/MooneyRivlin.h
  Integrators/InternalForce.h
  Integrators/MaterialTangent.h
  Integrators/BodyForce.h
  Integrators/TractionForce.h
  Constitutive/Hooke.h
  Fields/GreenLagrangeStrain.h
  Fields/FirstPiolaKirchhoffStress.h
  Fields/CauchyStress.h
  Linear/ForwardDecls.h
  Linear/LinearElasticityIntegral.h
  Linear/P1/LinearElasticityIntegral.h
)

# Header-only / interface library
add_library(RodinSolid INTERFACE)
add_library(Rodin::Solid ALIAS RodinSolid)
set_target_properties(RodinSolid PROPERTIES EXPORT_NAME Solid)

# ---- Includes / defs -------------------------------------------------------
target_include_directories(RodinSolid
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinSolid
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# ---- Link targets ----------------------------------------------------------
target_link_libraries(RodinSolid
  INTERFACE
    Rodin::Geometry
    Rodin::Variational
)

# Expose headers for IDEs in build tree
foreach(_h IN LISTS RodinSolid_HEADERS)
  target_sources(RodinSolid INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

# ---- Install target --------------------------------------------------------
install(
  TARGETS RodinSolid
  EXPORT RodinTargets
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Solid
)
