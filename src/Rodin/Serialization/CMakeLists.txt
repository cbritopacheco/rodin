set(RodinSerialization_HEADERS EigenMatrixSerialization.h)

set(RodinSerialization_SRCS )
add_library(RodinSerialization INTERFACE)
add_library(Rodin::Serialization ALIAS RodinSerialization)

target_link_libraries(RodinSerialization INTERFACE Boost::serialization)

