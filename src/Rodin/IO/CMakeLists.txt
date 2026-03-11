#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinIO_HEADERS
  MFEM.h
  MEDIT.h
  HDF5.h
  XDMF.h
  Loader.h
  Printer.h
  MeshLoader.h
  MeshPrinter.h
  GridFunctionLoader.h
  GridFunctionPrinter.h
)

set(RodinIO_SRCS
  MFEM.cpp
  MEDIT.cpp
  MeshLoader.cpp
  MeshPrinter.cpp
  GridFunctionLoader.cpp
)

add_library(RodinIO
  ${RodinIO_SRCS}
  ${RodinIO_HEADERS}
)
add_library(Rodin::IO ALIAS RodinIO)
set_target_properties(RodinIO PROPERTIES EXPORT_NAME IO)

# ---- Include dirs / defs ---------------------------------------------------
target_include_directories(RodinIO
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinIO
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# ---- Link targets ----------------------------------------------------------
target_link_libraries(RodinIO
  PUBLIC
    Rodin::Alert
    Rodin::Geometry
    Rodin::Variational
    Boost::filesystem
)

find_package(HDF5 COMPONENTS C QUIET)
if (HDF5_FOUND)
  target_compile_definitions(RodinIO PUBLIC RODIN_IO_HAS_HDF5=1)
  if (TARGET HDF5::HDF5)
    target_link_libraries(RodinIO PUBLIC HDF5::HDF5)
  else()
    target_include_directories(RodinIO PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(RodinIO PUBLIC ${HDF5_LIBRARIES})
  endif()
else()
  target_compile_definitions(RodinIO PUBLIC RODIN_IO_HAS_HDF5=0)
endif()

# ---- Install target --------------------------------------------------------
install(
  TARGETS RodinIO
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/IO
)
