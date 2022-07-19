add_executable(Trimming Trimming.cpp)
target_link_libraries(Trimming
  PUBLIC
  Rodin::Mesh)

add_executable(CCL CCL.cpp)
target_link_libraries(CCL
  PUBLIC
  Rodin::Mesh)

