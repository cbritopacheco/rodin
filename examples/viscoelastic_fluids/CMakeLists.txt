# Viscoelastic (sPTT/Oldroyd-B) flows through idealised arterial lesions.
# Same log-conformation formulation as Heart/LeftAtrium2D_viscous_logarithmic_implicit.
if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_executable(ArterialLesion2D_viscous_logarithmic_implicit
    ArterialLesion2D_viscous_logarithmic_implicit.cpp)
  target_link_libraries(ArterialLesion2D_viscous_logarithmic_implicit
    PUBLIC
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
  # OrthogonalProjection.h and CoronaryArtery/{Alerts,Timing}.h are shared with
  # the Heart examples.
  target_include_directories(ArterialLesion2D_viscous_logarithmic_implicit
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../Heart)
  if (RODIN_USE_SCOTCH)
    target_link_libraries(ArterialLesion2D_viscous_logarithmic_implicit PUBLIC Rodin::Scotch)
  endif()
endif()
