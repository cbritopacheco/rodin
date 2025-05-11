set(RodinIO_HEADERS
  MFEM.h
  MEDIT.h
  EnSight6.h
  Loader.h
  Printer.h
  MeshLoader.h
  MeshPrinter.h
  GridFunctionLoader.h
  GridFunctionPrinter.h)

set(RodinIO_SRCS
  MFEM.cpp
  MEDIT.cpp
  EnSight6.cpp
  MeshLoader.cpp
  MeshPrinter.cpp
  GridFunctionLoader.cpp)
add_library(RodinIO ${RodinIO_SRCS} ${RodinIO_HEADERS})
add_library(Rodin::IO ALIAS RodinIO)

target_link_libraries(RodinIO
  PUBLIC
  Rodin::Alert
  Rodin::Geometry
  Rodin::Variational
  Boost::filesystem)

install(TARGETS RodinIO
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/IO
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/IO
)
