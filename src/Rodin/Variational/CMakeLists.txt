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
  Problem.hpp
  ScalarCoefficient.hpp
  ScalarCoefficient.cpp
  NeumannBC.cpp
  DirichletBC.cpp
  VectorCoefficient.hpp
  VectorCoefficient.cpp
  DomainLFIntegrator.cpp
  DiffusionIntegrator.cpp
  ElasticityIntegrator.cpp
  )

add_library(RodinVariational
  ${RodinVariational_SRCS} ${RodinVariational_HEADERS})
add_library(Rodin::Variational ALIAS RodinVariational)

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinVariational
  INTERFACE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
  PUBLIC
  ${MFEM_INCLUDE_DIRS}
  )

target_link_libraries(RodinVariational
  PUBLIC
  RodinMesh
  RodinAlert
  RodinVariationalFormLanguage
  )

add_subdirectory(FormLanguage)
