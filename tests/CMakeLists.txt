find_package(Corrade REQUIRED TestSuite)
set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

enable_testing()

corrade_add_resource(
  RodinTest_RESOURCES "${PROJECT_SOURCE_DIR}/resources/Manifest.conf")

rodin_add_test(
  NAME Poisson
  MODULES Variational Geometry
  SOURCES Variational/Poisson.cpp ${RodinTest_RESOURCES}
  LIBRARIES Rodin::Geometry Rodin::Variational)
