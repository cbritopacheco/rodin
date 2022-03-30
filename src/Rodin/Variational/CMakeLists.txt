set(RodinVariational_HEADERS
  BilinearForm.h
  BilinearFormIntegrator.h
  DirichletBC.h
  EssentialBoundary.h
  Derivative.h
  Dot.h
  FiniteElementSpace.h
  FormLanguage.h
  ForwardDecls.h
  Grad.h
  GridFunction.h
  ShapeFunction.h
  TrialFunction.h
  TestFunction.h
  GridFunctionIndex.h
  GridFunctionView.h
  H1.h
  IdentityMatrix.h
  Integral.h
  Sum.h
  UnaryMinus.h
  Mult.h
  Jacobian.h
  LinearForm.h
  LinearFormIntegrator.h
  MatrixCoefficient.h
  Problem.h
  ProblemBody.h
  ScalarCoefficient.h
  Trace.h
  Transpose.h
  Rank3Operator.h
  Restriction.h
  VectorCoefficient.h
  )

set(RodinVariational_SRCS
  EssentialBoundary.cpp
  LinearFormIntegrator.cpp
  BilinearFormIntegrator.cpp
  ShapeFunction.cpp
  Derivative.cpp
  GridFunctionIndex.cpp
  GridFunctionView.cpp
  Integral.cpp
  Jacobian.cpp
  MatrixCoefficient.cpp
  Trace.cpp
  VectorCoefficient.cpp
  ScalarCoefficient.cpp
  Restriction.cpp
  Rank3Operator.cpp
  Sum.cpp
  UnaryMinus.cpp
  Mult.cpp
  BilinearForm.hpp
  Dot.cpp
  LinearForm.hpp
  Problem.hpp
  ProblemBody.cpp
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
  ${Boost_INCLUDE_DIRS}
  )

target_link_libraries(RodinVariational
  PUBLIC
  RodinMesh
  RodinAlert
  RodinUtility
  RodinVariationalFormLanguage
  )

add_subdirectory(FormLanguage)
