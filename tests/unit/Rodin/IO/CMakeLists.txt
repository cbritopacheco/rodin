add_executable(RodinIOMeshLoaderTest MeshLoaderTest.cpp)
target_link_libraries(RodinIOMeshLoaderTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinIOMeshLoaderTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinIOFileFormatTest FileFormatTest.cpp)
target_link_libraries(RodinIOFileFormatTest PRIVATE
  GTest::gtest_main Rodin::Geometry Rodin::Variational Rodin::IO)
gtest_discover_tests(RodinIOFileFormatTest PROPERTIES LABELS unit)

add_executable(RodinIOMEDITParserTest MEDITParserTest.cpp)
target_link_libraries(RodinIOMEDITParserTest PRIVATE GTest::gtest_main Rodin::IO)
gtest_discover_tests(RodinIOMEDITParserTest PROPERTIES LABELS unit)

add_executable(RodinIOH1GridFunctionIOTest H1GridFunctionIOTest.cpp)
target_link_libraries(RodinIOH1GridFunctionIOTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::IO)
gtest_discover_tests(RodinIOH1GridFunctionIOTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinIOP1GridFunctionIOTest P1GridFunctionIOTest.cpp)
target_link_libraries(RodinIOP1GridFunctionIOTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::IO)
gtest_discover_tests(RodinIOP1GridFunctionIOTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinIOP0GridFunctionIOTest P0GridFunctionIOTest.cpp)
target_link_libraries(RodinIOP0GridFunctionIOTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::IO)
gtest_discover_tests(RodinIOP0GridFunctionIOTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinIOHDF5IOTest HDF5IOTest.cpp)
target_link_libraries(RodinIOHDF5IOTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::IO)
gtest_discover_tests(RodinIOHDF5IOTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinIOSpatialMatrixTest SpatialMatrixTest.cpp)
target_link_libraries(RodinIOSpatialMatrixTest
  PUBLIC GTest::gtest GTest::gtest_main Rodin::IO Rodin::Variational)
gtest_discover_tests(RodinIOSpatialMatrixTest PROPERTIES LABELS unit)

add_executable(RodinIOXDMFTest XDMFTest.cpp)
target_link_libraries(RodinIOXDMFTest
  PUBLIC GTest::gtest GTest::gtest_main Rodin::IO Rodin::Variational)
gtest_discover_tests(RodinIOXDMFTest PROPERTIES LABELS unit)

add_executable(RodinIOMEDITLoaderTest MEDITLoaderTest.cpp)
target_link_libraries(RodinIOMEDITLoaderTest PRIVATE
  GTest::gtest_main Rodin::Geometry Rodin::Variational Rodin::IO)
gtest_discover_tests(RodinIOMEDITLoaderTest PROPERTIES LABELS unit)

add_executable(RodinIOMFEMCoefficientTest MFEMCoefficientTest.cpp)
target_link_libraries(RodinIOMFEMCoefficientTest PRIVATE
  GTest::gtest_main Rodin::Geometry Rodin::Variational Rodin::IO)
gtest_discover_tests(RodinIOMFEMCoefficientTest PROPERTIES LABELS unit)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_executable(RodinIOHDF5FailureTest HDF5FailureTest.cpp)
  target_link_libraries(RodinIOHDF5FailureTest PRIVATE
    GTest::gtest_main Rodin::Geometry Rodin::Variational Rodin::IO)
  foreach(operation IN ITEMS H5Screate_simple H5Dget_space
      H5Sget_simple_extent_ndims H5Sget_simple_extent_dims H5Dwrite H5Gcreate2)
    target_link_options(RodinIOHDF5FailureTest PRIVATE "LINKER:--wrap=${operation}")
  endforeach()
  gtest_discover_tests(RodinIOHDF5FailureTest PROPERTIES LABELS unit)
endif()
