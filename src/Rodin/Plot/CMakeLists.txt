# TODO: CURRENTLY BROKEN
set(RodinPlot_HEADERS )
set(RodinPlot_SRCS )

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

