#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]
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

# target_sources() does not resolve relative paths before 3.13
if (${CMAKE_VERSION} VERSION_LESS 3.13)
  list(
    TRANSFORM RodinVariationalFormLanguage_HEADERS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
endif()

add_library(RodinVariationalFormLanguage INTERFACE)
add_library(Rodin::Variational::FormLanguage ALIAS RodinVariationalFormLanguage)
target_sources(RodinVariationalFormLanguage
  INTERFACE ${RodinVariationalFormLanguage_HEADERS})

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinVariationalFormLanguage
  INTERFACE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
  )



