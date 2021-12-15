set(RodinVariational_HEADERS
  BilinearForm.h
  Derivative.h
  DiffusionIntegrator.h
  DirichletBC.h
  DomainLFIntegrator.h
  ElasticityIntegrator.h
  FiniteElementSpace.h
  FormLanguage.h
  ForwardDecls.h
  Jacobian.h
  GridFunction.h
  H1.h
  L2.h
  LinearForm.h
  MatrixCoefficient.h
  NeumannBC.h
  Problem.h
  ScalarCoefficient.h
  Transpose.h
  VectorCoefficient.h
  Dot.h
  Trace.h
  )

set(RodinVariational_SRCS
  Derivative.hpp
  Problem.hpp
  ScalarCoefficient.hpp
  VectorCoefficient.hpp
  Dot.hpp
  DiffusionIntegrator.cpp
  DirichletBC.cpp
  DomainLFIntegrator.cpp
  ElasticityIntegrator.cpp
  NeumannBC.cpp
  Transpose.cpp
  VectorCoefficient.cpp
  MatrixCoefficient.cpp
  Jacobian.cpp
  Trace.cpp
  Gradient.cpp
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
