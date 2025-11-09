set(RodinPlot_HEADERS
  Artist/Axes/Axes2D.h
  Artist/Axes/Axes3D.h
  Artist/Axis/Axis.h
  Artist/Figure.h
  Artist/Lines/Line2D.h
  Backend/Bases/BaseArtist2D.h
  Backend/Bases/BaseAxes.h
  Backend/Bases/BaseCamera2D.h
  Backend/Bases/BaseCamera3D.h
  Backend/Bases/BaseDrawable2D.h
  Backend/Bases/BaseDrawable3D.h
  Backend/Bases/BaseEvent.h
  Backend/Bases/BaseFigure.h
  Backend/Bases/BaseTopLevelArtist2D.h
  Backend/Bases/FigureCanvas.h
  Backend/Event/MouseEvent/MouseButtonEvent.h
  Backend/Event/MouseEvent/MouseEvent.h
  Backend/Event/MouseEvent/MouseMotionEvent.h
  Backend/Event/MouseEvent/MouseWheelEvent.h
  Backend/Event/MouseEvent.h
  Backend/Event.h
  Backend/Renderer/Camera2D.h
  Backend/Renderer/Common.h
  Backend/Renderer/Drawables/Frame.h
  Backend/Renderer/Drawables/Line2D.h
  Common.h
  Configuration.h
  ForwardDecls.h
  GUI/Cursor.h
  Plot.h
  )

set(RodinPlot_SRCS
  Artist/Axes/Axes2D.hpp
  Artist/Figure.hpp
  Artist/Axes/Axes2D.cpp
  Artist/Axis/Axis.cpp
  Artist/Figure.cpp
  Backend/Bases/BaseArtist2D.hpp
  Backend/Bases/BaseArtist2D.cpp
  Backend/Bases/BaseAxes.cpp
  Backend/Bases/BaseFigure.cpp
  Backend/Bases/BaseTopLevelArtist2D.cpp
  Backend/Renderer/Drawables/Line2D.cpp
  Configuration.cpp
  Plot.cpp
  )

add_library(RodinPlot ${RodinPlot_SRCS} ${RodinPlot_HEADERS} )
add_library(Rodin::Plot ALIAS RodinPlot)
set_target_properties(RodinPlot PROPERTIES EXPORT_NAME Plot)

target_include_directories(RodinPlot
  INTERFACE
  $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
  )

# ---- Link targets ----------------------------------------------------------
target_link_libraries(RodinPlot
  PUBLIC
  Rodin::Rodin
  Eigen3::Eigen
  Corrade::Main
  Magnum::GL
  Magnum::Magnum
  Magnum::Shaders
  Magnum::Primitives
  Magnum::SceneGraph
  Magnum::Application
  MagnumIntegration::Eigen
  )

install(TARGETS RodinPlot
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/Plot
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Plot
)
