add_executable(AnnulusDetection AnnulusDetection.cpp)
target_link_libraries(AnnulusDetection PUBLIC Rodin::Geometry)

add_executable(AnnulusClosing AnnulusClosing.cpp)
target_link_libraries(AnnulusClosing PUBLIC Rodin::Geometry)

add_executable(AttributeDetection AttributeDetection.cpp)
target_link_libraries(AttributeDetection PUBLIC Rodin::Geometry)
