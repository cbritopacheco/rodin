#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinT8Code_HEADERS
  Mesh.h
  ForwardDecls.h
)

set(RodinT8Code_SRCS
  Mesh.cpp
)

add_library(RodinT8Code
  ${RodinT8Code_HEADERS} ${RodinT8Code_SRCS})
add_library(Rodin::T8Code ALIAS RodinT8Code)

target_link_libraries(RodinT8Code
  PUBLIC
  T8Code::T8Code
  Rodin::Alert
  Rodin::Geometry
)

# Note: RodinT8Code is not exported because it depends on t8code which is an
# external library. Users who want t8code support should build Rodin from
# source with -DRODIN_USE_T8CODE=ON.
# install(TARGETS RodinT8Code
#   EXPORT       RodinTargets
#   RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
#   LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
#   ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
#   INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/T8Code
# )
