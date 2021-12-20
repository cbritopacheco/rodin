set(RodinVariationalFormLanguage_HEADERS
  Base.h
  BilinearFormIntegratorSum.h
  BilinearFormIntegratorUnaryMinus.h
  BoundaryConditionList.h
  ForwardDecls.h
  LinearFormIntegratorSum.h
  LinearFormIntegratorUnaryMinus.h
  MatrixSum.h
  ProblemBody.h
  ScalarMatrixProduct.h
  ScalarProduct.h
  ScalarSum.h
  ScalarUnaryMinus.h
  )

set(RodinVariationalFormLanguage_SRCS
  BilinearFormIntegratorSum.cpp
  BilinearFormIntegratorUnaryMinus.cpp
  BoundaryConditionList.cpp
  LinearFormIntegratorSum.cpp
  LinearFormIntegratorUnaryMinus.cpp
  MatrixSum.cpp
  ProblemBody.cpp
  ScalarMatrixProduct.cpp
  ScalarProduct.cpp
  ScalarSum.cpp
  ScalarUnaryMinus.cpp
  )

add_library(RodinVariationalFormLanguage
  ${RodinVariationalFormLanguage_SRCS} ${RodinVariationalFormLanguage_HEADERS})
add_library(Rodin::Variational::FormLanguage ALIAS RodinVariationalFormLanguage)

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinVariationalFormLanguage
  PUBLIC
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
  PRIVATE
  ${MFEM_INCLUDE_DIRS}
  )



