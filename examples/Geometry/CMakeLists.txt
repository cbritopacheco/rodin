add_executable(Trace Trace.cpp)
target_link_libraries(Trace
  PUBLIC
  Rodin::Geometry)

add_executable(Volume Volume.cpp)
target_link_libraries(Volume
  PUBLIC
  Rodin::Geometry)

add_executable(Displace Displace.cpp)
target_link_libraries(Displace
  PUBLIC
  Rodin::Geometry
  Rodin::Variational)

add_executable(Skin Skin.cpp)
target_link_libraries(Skin
  PUBLIC
  Rodin::Geometry)

add_executable(Trim Trim.cpp)
target_link_libraries(Trim
  PUBLIC
  Rodin::Geometry)

add_executable(CCL CCL.cpp)
target_link_libraries(CCL
  PUBLIC
  Rodin::Geometry)

add_executable(Connectivity Connectivity.cpp)
target_link_libraries(Connectivity
  PUBLIC
  Rodin::Geometry)

add_executable(Square Square.cpp)
target_link_libraries(Square
  PUBLIC
  Rodin::Geometry)

add_executable(UniformGrid UniformGrid.cpp)
target_link_libraries(UniformGrid
  PUBLIC
  Rodin::Geometry)

add_executable(Box Box.cpp)
target_link_libraries(Box
  PUBLIC
  Rodin::Geometry)

add_executable(GreedyPartitioner GreedyPartitioner.cpp)
target_link_libraries(GreedyPartitioner
  PUBLIC
  Rodin::Geometry)

add_executable(BalancedCompactPartitioner BalancedCompactPartitioner.cpp)
target_link_libraries(BalancedCompactPartitioner
  PUBLIC
  Rodin::Geometry)

add_executable(Scale Scale.cpp)
target_link_libraries(Scale
  PUBLIC
  Rodin::Geometry)

add_executable(AnnulusDetection AnnulusDetection.cpp)
target_link_libraries(AnnulusDetection PUBLIC Rodin::Geometry)

add_executable(AnnulusClosing AnnulusClosing.cpp)
target_link_libraries(AnnulusClosing PUBLIC Rodin::Geometry)

add_executable(AttributeDetection AttributeDetection.cpp)
target_link_libraries(AttributeDetection PUBLIC Rodin::Geometry)
