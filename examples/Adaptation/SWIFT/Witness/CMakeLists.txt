option(RODIN_BUILD_SWIFT_WITNESS "Build the offline reference witness generator" OFF)

if(RODIN_BUILD_SWIFT_WITNESS)
  find_package(Boost 1.75 REQUIRED COMPONENTS json)
  find_package(CGAL REQUIRED)
  add_executable(SWIFT_Witness Witness.cpp)
  target_link_libraries(SWIFT_Witness PRIVATE Rodin::Geometry Boost::json CGAL::CGAL)
endif()
