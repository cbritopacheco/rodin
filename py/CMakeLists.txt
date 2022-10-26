set(PyRodin_SRCS rodin.cpp)

add_library(PyRodin MODULE ${PyRodin_SRCS})
set_property(TARGET PyRodin PROPERTY OUTPUT_NAME rodin)
target_link_libraries(PyRodin
  PUBLIC
  Rodin::Geometry
  PRIVATE
  pybind11::module)
pybind11_extension(PyRodin)
target_compile_definitions(PyRodin PRIVATE VERSION_INFO=${PROJECT_VERSION})

