set(RodinVariationalFormLanguage_HEADERS
  BCExpr.h
  BCExprList.h
  BCExprListCons.h
  BFExprBCExprSum.h
  BFExprLFExprSum.h
  BilinearFormExpr.h
  BilinearFormExprSum.h
  BilinearFormExprUnaryMinus.h
  CMakeLists.txt
  FormExprBCExprSum.h
  ForwardDecls.h
  LinearFormExpr.h
  LinearFormExprSum.h
  LinearFormExprUnaryMinus.h
  ProblemBody.h
  RodinBase.h
  ScalarCoefficientSum.h
  ScalarCoefficientUnaryMinus.h
  TypeTraits.h
  )

set(RodinVariationalFormLanguage_SRCS
  ScalarCoefficientUnaryMinus.cpp
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



