#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

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
  IdentityMatrix.h
  Integral.h
  Sum.h
  UnaryMinus.h
  Mult.h
  Jacobian.h
  LinearForm.h
  LinearFormIntegrator.h
  MatrixFunction.h
  Problem.h
  ProblemBody.h
  RealFunction.h
  Trace.h
  Transpose.h
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
  LinearElasticity/LinearElasticityIntegral.h
)

set(RodinVariational_SRCS
  P0/P0.cpp
  P0/P0Element.cpp
  P1/P1.cpp
  P1/P1Element.cpp
  H1/H1.cpp
  H1/H1Element.cpp
  )

add_library(RodinVariational
  ${RodinVariational_SRCS}
  ${RodinVariational_HEADERS}
)
add_library(Rodin::Variational ALIAS RodinVariational)
set_target_properties(RodinVariational PROPERTIES EXPORT_NAME Variational)

# ---- Include / compile interface -------------------------------------------
target_include_directories(RodinVariational
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinVariational
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# ---- Link targets ----------------------------------------------------------
target_link_libraries(RodinVariational
  PUBLIC
    Rodin::QF
    Rodin::IO
    Rodin::Math
    Rodin::Alert
    Rodin::Utility
    Rodin::Assembly
    Rodin::Geometry
    Rodin::FormLanguage
)

# ---- Install ---------------------------------------------------------------
install(TARGETS RodinVariational
  EXPORT  RodinTargets
  RUNTIME DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Variational
)
