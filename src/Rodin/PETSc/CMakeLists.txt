#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]
set(RodinPETSc_HEADERS
  ForwardDecls.h
  Object.h
  IO/GridFunctionPrinter.h
  Solver/KSP.h
  Solver/CG.h
  Math/Vector.h
  Math/Matrix.h
  Math/LinearSystem.h
  Assembly/MPI.h
  Assembly/Sequential.h
  Variational/GridFunction.h
  )

set(RodinPETSc_SRCS
  IO/GridFunctionPrinter.cpp
  Solver/KSP.cpp
  Solver/CG.cpp
  Assembly/MPI.cpp
  Assembly/OpenMP.cpp
  Assembly/Sequential.cpp
  Math/LinearSystem.cpp
  Variational/GridFunction.cpp
  Variational/TrialFunction.cpp
  )

add_library(RodinPETSc ${RodinPETSc_HEADERS} ${RodinPETSc_SRCS})
add_library(Rodin::PETSc ALIAS RodinPETSc)

target_include_directories(RodinPETSc
  INTERFACE $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>)

find_package(PETSc REQUIRED)

target_link_libraries(RodinPETSc
  PUBLIC
  Rodin::Geometry
  Rodin::Alert
  Rodin::MPI
  Rodin::Math
  Rodin::Variational
  Rodin::Assembly
  Rodin::FormLanguage
  PETSc::PETSc)

if (RODIN_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(RodinPETSc PUBLIC OpenMP::OpenMP_CXX)
endif()

