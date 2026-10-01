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

# The same formulation in 3D (pipe), in its own folder.
add_subdirectory(ArterialLesion3D_viscous_logarithmic_implicit)

# Convergence tests of the same solver (Oldroyd-B), meshes built in the code.
if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_executable(test_viscoelastic_implicit test_viscoelastic_implicit.cpp)
  target_link_libraries(test_viscoelastic_implicit
    PUBLIC
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
  target_include_directories(test_viscoelastic_implicit
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../Heart)
  if (RODIN_USE_SCOTCH)
    target_link_libraries(test_viscoelastic_implicit PUBLIC Rodin::Scotch)
  endif()
endif()

# Viscoelastic Kolmogorov flow (Berti & Boffetta 2010), doubly periodic box.
if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_executable(test_kolmogorov_viscoelastic test_kolmogorov_viscoelastic.cpp)
  target_link_libraries(test_kolmogorov_viscoelastic
    PUBLIC
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
  target_include_directories(test_kolmogorov_viscoelastic PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endif()
