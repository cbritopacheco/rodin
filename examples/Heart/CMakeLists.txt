add_executable(CCMLC2014_0D CCMLC2014_0D.cpp)
target_link_libraries(CCMLC2014_0D
  PUBLIC
  Rodin::Heart
  Rodin::Alert
  )

add_executable(PoroelasticSphere_0D PoroelasticSphere_0D.cpp)
target_link_libraries(PoroelasticSphere_0D
  PUBLIC
  Rodin::Heart
  Rodin::Alert
  )

add_executable(CoronarySolid CoronarySolid.cpp)
target_link_libraries(CoronarySolid
  PUBLIC
  Rodin::Solid
  Rodin::Heart
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  )

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_executable(CoronaryArterySolid_PETSc_MPI CoronaryArterySolid_PETSc_MPI.cpp)
  target_link_libraries(CoronaryArterySolid_PETSc_MPI
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
Rodin::MPI
    Rodin::Alert Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
  if (RODIN_USE_SCOTCH)
    target_link_libraries(CoronaryArterySolid_PETSc_MPI PUBLIC Rodin::Scotch)
  endif()

  add_executable(CoronaryArtery_FSI_PETSc_MPI CoronaryArtery_FSI_PETSc_MPI.cpp)
  target_link_libraries(CoronaryArtery_FSI_PETSc_MPI
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
  if (RODIN_USE_SCOTCH)
    target_link_libraries(CoronaryArtery_FSI_PETSc_MPI PUBLIC Rodin::Scotch)
  endif()

  add_executable(CoronaryArtery_FSI_Explicit_PETSc_MPI CoronaryArtery_FSI_Explicit_PETSc_MPI.cpp)
  target_link_libraries(CoronaryArtery_FSI_Explicit_PETSc_MPI
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
  if (RODIN_USE_SCOTCH)
    target_link_libraries(CoronaryArtery_FSI_Explicit_PETSc_MPI PUBLIC Rodin::Scotch)
  endif()
endif()

# Sequential (single-process) explicit coronary FSI: PETSc only, no MPI.
if (RODIN_USE_PETSC)
  add_executable(CoronaryArtery_FSI_Explicit_PETSc_Seq CoronaryArtery_FSI_Explicit_PETSc_Seq.cpp)
  target_link_libraries(CoronaryArtery_FSI_Explicit_PETSc_Seq
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)

  # Equal-order P1/P1 variant (separate results_p1p1/ output folder), so it
  # can run concurrently with the P2/P1 case on different cores.
  add_executable(CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq.cpp)
  target_link_libraries(CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)

  # Closed-loop coronary FSI: P1/P1 explicit ALE FSI coupled (staggered) to
  # the 0D poroelastic left ventricle.  Sequential, PETSc only.
  add_executable(CoronaryArtery_Explicit_PoroElastic CoronaryArtery_Explicit_PoroElastic.cpp)
  target_link_libraries(CoronaryArtery_Explicit_PoroElastic
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)

  # Coagulation-kinetics testbed on a simplified geometry: one inlet, one
  # outlet, an R_mu-Windkessel outlet and an ectatic segment that provides the
  # recirculation the kinetics need. Sequential and PETSc-only, like the P1/P1
  # case above; it does not link Rodin::MPI.
  #
  # The mesh is not in the repository because it is generated:
  #   python3 examples/Heart/ProteinTest/make_test_mesh.py \
  #           examples/Heart/ProteinTest/ProteinTest.mesh
  # and the executable looks for it at ProteinTest/ProteinTest.mesh relative to
  # the working directory, so run it from the examples/Heart tree or pass a
  # different path in Config::meshPath.
  add_executable(ProteinTest ProteinTest.cpp)
  target_link_libraries(ProteinTest
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
  target_include_directories(ProteinTest
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endif()

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  # Atrium: four pressure inlets, one outlet closed by the 0D LV through a
  # mitral diode, VMS/PSPG-stabilised P1/P1 flow and the coagulation kinetics.
  add_executable(Atrium Atrium.cpp)
  target_link_libraries(Atrium
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
  target_include_directories(Atrium PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  if (RODIN_USE_SCOTCH)
    target_link_libraries(Atrium PUBLIC Rodin::Scotch)
  endif()

  # Two-dimensional left atrium with a rectangular appendage. Same stabilised
  # P1/P1 formulation and coagulation kinetics as Atrium, but driven by the two
  # measured pressure waveforms in resources/examples/Heart rather than by a 0D
  # ventricle. The mesh must be flattened once before the first run:
  #   python3 examples/Heart/LA2D/make_la2d_mesh.py \
  #           resources/examples/Heart/LA2D_rectLAA.mesh \
  #           resources/examples/Heart/LA2D_rectLAA_2D.mesh
  # LA2D_rectLAA.mesh is stored as MEDIT "Dimension 3" with z = 0, and Rodin
  # reads the space dimension from that keyword, so loading it unchanged would
  # give a three-component velocity on a surface.
  add_executable(LeftAtrium2D LeftAtrium2D.cpp)
  target_link_libraries(LeftAtrium2D
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
  target_include_directories(LeftAtrium2D PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  if (RODIN_USE_SCOTCH)
    target_link_libraries(LeftAtrium2D PUBLIC Rodin::Scotch)
  endif()

  add_executable(CoronaryArtery
    CoronaryArtery.cpp
    CoronaryArtery/CoupledLV0DCoronary3D.cpp
    )
  target_link_libraries(CoronaryArtery
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    )

  if (RODIN_USE_SCOTCH)
    target_link_libraries(CoronaryArtery PUBLIC Rodin::Scotch)
  endif()
endif()
