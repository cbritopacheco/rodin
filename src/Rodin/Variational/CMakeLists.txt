set(RodinVariational_HEADERS
  Advect.h
  BilinearForm.h
  BilinearFormIntegrator.h
  BoundaryCondition.h
  Derivative.h
  DiffusionIntegrator.h
  DirichletBC.h
  DomainLFIntegrator.h
  Dot.h
  ElasticityIntegrator.h
  FiniteElementSpace.h
  FormLanguage.h
  ForwardDecls.h
  Gradient.h
  GridFunction.h
  H1.h
  IdentityMatrix.h
  Jacobian.h
  LinearForm.h
  LinearFormIntegrator.h
  MatrixCoefficient.h
  NeumannBC.h
  Problem.h
  ScalarCoefficient.h
  Trace.h
  Transpose.h
  VectorBoundaryFluxLFIntegrator.h
  VectorCoefficient.h
  MassIntegrator.h
  VectorMassIntegrator.h
  VectorDiffusionIntegrator.h
  VectorDomainLFIntegrator.h
  )

set(RodinVariational_SRCS
  Advect.cpp
  DiffusionIntegrator.cpp
  DirichletBC.cpp
  DomainLFIntegrator.cpp
  ElasticityIntegrator.cpp
  Gradient.cpp
  IdentityMatrix.cpp
  Jacobian.cpp
  MatrixCoefficient.cpp
  NeumannBC.cpp
  Trace.cpp
  Transpose.cpp
  VectorBoundaryFluxLFIntegrator.cpp
  VectorCoefficient.cpp
  MassIntegrator.cpp
  VectorMassIntegrator.cpp
  VectorDiffusionIntegrator.cpp
  VectorDomainLFIntegrator.cpp
  ScalarCoefficient.cpp
  BilinearForm.hpp
  Derivative.hpp
  Dot.hpp
  LinearForm.hpp
  Problem.hpp
  ScalarCoefficient.hpp
  VectorCoefficient.hpp
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
