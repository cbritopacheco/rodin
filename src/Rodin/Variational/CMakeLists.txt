set(RodinVariational_HEADERS
  FiniteElement.h
  BilinearForm.h
  BilinearFormIntegrator.h
  DirichletBC.h
  Derivative.h
  Dot.h
  Jump.h
  FiniteElementSpace.h
  ForwardDecls.h
  Grad.h
  Function.h
  GridFunction.h
  ShapeFunction.h
  TrialFunction.h
  TestFunction.h
  H1.h
  IdentityMatrix.h
  Integral.h
  Sum.h
  MFEM.h
  UnaryMinus.h
  Mult.h
  Jacobian.h
  LinearForm.h
  LinearFormIntegrator.h
  MatrixFunction.h
  Problem.h
  ProblemBody.h
  ScalarFunction.h
  Trace.h
  Transpose.h
  TensorBasis.h
  Restriction.h
  VectorFunction.h
  Minus.h
  Min.h
  Max.h
  GEQ.h
  LEQ.h
  GT.h
  LT.h
  OR.h
  AND.h
  Tangent.h
  Assembly/AssemblyBase.h
  Assembly/Native.h
  LinearElasticity/LinearElasticityIntegral.h
  )

set(RodinVariational_SRCS
  FiniteElement.cpp
  FiniteElementSpace.cpp
  LinearFormIntegrator.cpp
  BilinearFormIntegrator.cpp
  ShapeFunction.cpp
  Derivative.cpp
  Integral.cpp
  Jacobian.cpp
  GridFunction.cpp
  Division.cpp
  Grad.cpp
  MatrixFunction.cpp
  Trace.cpp
  VectorFunction.cpp
  ScalarFunction.cpp
  Restriction.cpp
  Sum.cpp
  UnaryMinus.cpp
  Mult.cpp
  BilinearForm.hpp
  Dot.cpp
  LinearForm.hpp
  Problem.hpp
  ProblemBody.cpp
  MFEM.cpp
  Minus.cpp
  Function.cpp
  RangeShape.cpp
  RangeType.cpp
  Composition.cpp
  Max.cpp
  Min.cpp
  GEQ.cpp
  LEQ.cpp
  GT.cpp
  LT.cpp
  OR.cpp
  QuadratureRule.cpp
  AND.cpp
  Tangent.cpp
  Assembly/AssemblyBase.cpp
  Assembly/Native.cpp
  )

add_library(RodinVariational
  ${RodinVariational_SRCS} ${RodinVariational_HEADERS})
add_library(Rodin::Variational ALIAS RodinVariational)

if (RODIN_USE_OPENMP)
  target_sources(RodinVariational PUBLIC Assembly/OpenMP.h Assembly/OpenMP.cpp)
  target_link_libraries(RodinVariational PUBLIC OpenMP::OpenMP_CXX)
endif()

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinVariational
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

target_link_libraries(RodinVariational
  PRIVATE
  mfem
  PUBLIC
  Rodin::IO
  Rodin::Alert
  Rodin::Utility
  Rodin::Geometry
  Rodin::FormLanguage)
