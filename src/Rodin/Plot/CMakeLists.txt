#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]
set(RodinPlot_HEADERS "")
set(RodinPlot_SRCS "")

add_library(RodinPlot ${RodinPlot_HEADERS} ${RodinPlot_SRCS})
add_library(Rodin::Plot ALIAS RodinPlot)

# ---- Link targets ----------------------------------------------------------
target_link_libraries(Rodin::Plot
  PUBLIC
  Eigen3::Eigen
  Corrade::Main
  Magnum::GL
  Magnum::Magnum
  Magnum::Shaders
  Magnum::Primitives
  Magnum::SceneGraph
  Magnum::Application
  PRIVATE
  MagnumIntegration::Eigen
  )

