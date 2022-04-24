set(RodinCoreGeometry_HEADERS
  Base.h
  Circle.h
  Concepts.h
  ForwardDecls.h
  Line2D.h
  LineSegment2D.h
  Point2D.h
  Rectangle.h)

set(RodinCoreGeometry_SRCS
  Circle.hpp
  Line2D.hpp
  LineSegment2D.hpp
  Point2D.hpp
  Rectangle.hpp)

# target_sources() does not resolve relative paths before 3.13
if (${CMAKE_VERSION} VERSION_LESS 3.13)
  list(
    TRANSFORM RodinCoreGeometry_SRCS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
  list(
    TRANSFORM RodinCoreGeometry_HEADERS
    PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/"
    )
endif()

add_library(RodinCoreGeometry INTERFACE)
add_library(Rodin::Core::Geometry ALIAS RodinCore)
target_sources(RodinCoreGeometry
  INTERFACE
  ${RodinCoreGeometry_SRCS}
  ${RodinCoreGeometry_HEADERS})

# ---- Link targets ----------------------------------------------------------
target_link_libraries(RodinCoreGeometry
  PUBLIC
  Rodin::Core
  Eigen3::Eigen)
