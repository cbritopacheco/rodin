set(RodinVariationalFormLanguage_HEADERS
  Base.h
  ForwardDecls.h
  )

set(RodinVariationalFormLanguage_SRCS
  Base.cpp
  )

add_library(RodinVariationalFormLanguage
  ${RodinVariationalFormLanguage_SRCS} ${RodinVariationalFormLanguage_HEADERS})
add_library(Rodin::FormLanguage ALIAS RodinVariationalFormLanguage)

# ---- Link targets ----------------------------------------------------------
target_include_directories(RodinVariationalFormLanguage
  PUBLIC
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
  PRIVATE
  ${MFEM_INCLUDE_DIRS}
  )

target_link_libraries(RodinVariationalFormLanguage
  PUBLIC
  RodinVariational
  )



