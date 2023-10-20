set(RodinMath_HEADERS
  Constants.h
  Common.h
  Vector.h
  Matrix.h)

set(RodinMath_SRCS
  Common.cpp
  Matrix.cpp
  SparseMatrix.cpp)

# target_sources() does not resolve relative paths before 3.13
if (${CMAKE_VERSION} VERSION_LESS 3.13)
  list(
    TRANSFORM RodinMath_SRCS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
  list(
    TRANSFORM RodinMath_HEADERS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
endif()

add_library(RodinMath ${RodinMath_SRCS} ${RodinMath_HEADERS})
add_library(Rodin::Math ALIAS RodinMath)

# ---- Link targets ----------------------------------------------------------
target_link_libraries(RodinMath
  PUBLIC
  Rodin
  Eigen3::Eigen)

target_include_directories(RodinMath
  INTERFACE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

