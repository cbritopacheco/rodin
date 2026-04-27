add_executable(CCMLC2014_0D CCMLC2014_0D.cpp)
target_link_libraries(CCMLC2014_0D
  PUBLIC
  Rodin::Heart
  Rodin::Alert
  )

add_executable(CoronaryArtery
  CoronaryArtery.cpp
  CoronaryArtery/CoupledLV0DCoronary3D.cpp
  )
target_link_libraries(CoronaryArtery
  PUBLIC
  Rodin::Heart
  Rodin::Variational
  Rodin::PETSc
  Rodin::Alert
  )
