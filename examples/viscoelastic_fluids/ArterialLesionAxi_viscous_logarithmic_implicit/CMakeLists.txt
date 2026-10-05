# Axisymmetric counterpart of ArterialLesion2D_viscous_logarithmic_implicit:
# the same log-conformation sPTT/Oldroyd-B formulation on the meridian
# half-plane of the pipe (hoop terms, r-weighted forms).
if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_executable(ArterialLesionAxi_viscous_logarithmic_implicit
    ArterialLesionAxi_viscous_logarithmic_implicit.cpp)
  target_link_libraries(ArterialLesionAxi_viscous_logarithmic_implicit
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
  target_include_directories(ArterialLesionAxi_viscous_logarithmic_implicit
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../Heart)
  if (RODIN_USE_SCOTCH)
    target_link_libraries(ArterialLesionAxi_viscous_logarithmic_implicit PUBLIC Rodin::Scotch)
  endif()
endif()
