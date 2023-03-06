set(RodinUtility_HEADERS
  Cast.h
  IsSpecialization.h
  OptionalReference.h
  Overloaded.h
  )

set(RodinUtility_SRCS Cast.cpp)

add_library(RodinUtility ${RodinUtility_SRCS} ${RodinUtility_HEADERS})
add_library(Rodin::Utility ALIAS RodinUtility)

target_include_directories(RodinUtility
  INTERFACE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
  PUBLIC
  ${MFEM_INCLUDE_DIRS})

