#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]
set(RodinVariational_HEADERS
  BilinearForm.h
  DiffusionIntegrator.h
  DirichletBC.h
  DomainLFIntegrator.h
  ElasticityIntegrator.h
  FiniteElementSpace.h
  FormLanguage.h
  ForwardDecls.h
  GridFunction.h
  H1.h
  LinearForm.h
  NeumannBC.h
  Problem.h
  ScalarCoefficient.h
  VectorCoefficient.h
  )

set(RodinVariational_SRCS
  DiffusionIntegrator.hpp
  DirichletBC.hpp
  DomainLFIntegrator.hpp
  ElasticityIntegrator.hpp
  NeumannBC.hpp
  Problem.hpp
  ScalarCoefficient.hpp
  VectorCoefficient.hpp
  )

# target_sources() does not resolve relative paths before 3.13
if (${CMAKE_VERSION} VERSION_LESS 3.13)
  list(
    TRANSFORM RodinVariational_SRCS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
  list(
    TRANSFORM RodinVariational_HEADERS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
endif()

add_library(RodinVariational INTERFACE)
target_sources(RodinVariational
  INTERFACE ${RodinVariational_SRCS} ${RodinVariational_HEADERS})

add_library(Rodin::Variational ALIAS RodinVariational)

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinVariational
  INTERFACE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
  ${MFEM_INCLUDE_DIRS}
  )

target_link_libraries(RodinVariational
  INTERFACE
  RodinMesh
  RodinAlert
  RodinVariationalFormLanguage
  )

add_subdirectory(FormLanguage)
