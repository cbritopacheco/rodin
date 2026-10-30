# TODO: migration du build Kokkos + PETSc sur demeter

Machine locale (Mint/Jammy) abandonnée: GCC 11 seulement, et Rodin exige
GCC 12+ ou Clang 14+ (`doc/Guides/Installation.dox`). La CI tourne sur
ubuntu-24.04 avec gcc-14. On passe tout sur demeter (gcc 14 natif, plus de
coeurs).

## Etat actuel

- PR #379 (fixes CUDA du module Kokkos, cible `module/Kokkos`): mergée.
- `tests/unit/Rodin/Math/SpatialComplexTest.cpp`: modif locale non commitée
  (`left({1, 1, 1})` -> `left(1, 1, 1)`), uniquement pour contourner GCC 11.
  A jeter sur demeter si le test compile tel quel avec gcc 14.
- Autres erreurs GCC 11 vues dans les tests, a revérifier avec gcc 14
  (normalement elles disparaissent):
  - `HingeProblem.h:28` et `:124` (RealFunction depuis une lambda),
    utilisé par SWIFTAdaptTest, SWIFTSolverTest, SWIFTExampleParametersTest;
  - `DirichletBCTest.cpp:212` et `:436` (CTAD sur RealFunction);
  - `IdentityTransformationTest.cpp:107-108` (destructeur virtuel constexpr
    dans un static_assert).

## 1. Vérifier l'environnement demeter

- [ ] `g++-14 --version`, `gfortran-14 --version` (gfortran requis pour MUMPS
      et ScaLAPACK).
- [ ] GPU et CUDA: `nvidia-smi`, `nvcc --version`. Noter l'arch
      (sm_XX) pour `Kokkos_ARCH_*` et `--with-cuda-arch`.
- [ ] MPI: quelle implémentation, quelle version, CUDA-aware ou non
      (`ompi_info --parsable --all | grep mpi_built_with_cuda_support`).
      Un seul MPI dans le PATH, sinon FindMPI de Rodin et PETSc divergent.
- [ ] CMake: si CMake 4, exporter `CMAKE_POLICY_VERSION_MINIMUM=3.5`
      (sinon METIS ne configure pas).
- [ ] `nproc` et RAM: choisir les `-j` (Rodin consomme beaucoup de RAM par
      unité de compilation, compter environ 3-4 Go par job).

## 2. Kokkos (build partagé, PIC)

PETSc en bibliothèque partagée exige un Kokkos partagé. Rodin doit utiliser
exactement le même Kokkos que PETSc (sinon deux runtimes Kokkos).

```sh
K=$HOME/kokkos-src   # cloner kokkos 4.7.x ici
export NVCC_WRAPPER_DEFAULT_COMPILER=g++-14
cmake -S $K/kokkos -B $K/build-shared \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=$K/kokkos/bin/nvcc_wrapper \
  -DCMAKE_CXX_STANDARD=20 \
  -DCMAKE_INSTALL_PREFIX=$K/install-shared \
  -DBUILD_SHARED_LIBS=ON \
  -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_SERIAL=ON \
  -DKokkos_ARCH_<ARCH>=ON \
  -DKokkos_ENABLE_IMPL_VIEW_LEGACY=ON \
  -DKokkos_ENABLE_COMPLEX_ALIGN=ON \
  -DKokkos_ENABLE_DEPRECATION_WARNINGS=ON
cmake --build $K/build-shared -j && cmake --install $K/build-shared
```

- [ ] Build et install de Kokkos partagé.
- [ ] Optionnel: même chose en statique (`build`/`install`) pour les builds
      Rodin sans PETSc.

## 3. PETSc 3.26 avec Kokkos, CUDA, hypre, SuperLU_DIST, MUMPS

PETSc 3.19 ne supporte pas Kokkos 4.7: prendre 3.26 avec Kokkos Kernels
4.7.04. PETSc n'a pas d'interface cuDSS (même en 3.26).

```sh
export OMPI_CC=gcc-14 OMPI_CXX=g++-14 OMPI_FC=gfortran-14   # si OpenMPI
./configure PETSC_ARCH=arch-cuda-kokkos \
  --prefix=$HOME/petsc-kokkos-install \
  --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
  --with-cxx-dialect=20 --with-debugging=0 \
  COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3 \
  --with-cuda=1 --with-cuda-dir=<CUDA_DIR> --with-cuda-arch=<SM> \
  --with-kokkos-dir=$K/install-shared \
  --download-kokkos-kernels --download-kokkos-kernels-commit=4.7.04 \
  --download-fblaslapack --download-scalapack --download-mumps \
  --download-metis --download-parmetis --download-superlu_dist \
  --download-hypre
make PETSC_DIR=$PWD PETSC_ARCH=arch-cuda-kokkos all
make PETSC_DIR=$PWD PETSC_ARCH=arch-cuda-kokkos install
```

- [ ] Configure + build + install.
- [ ] Si MPI n'est pas CUDA-aware: `export PETSC_OPTIONS=-use_gpu_aware_mpi 0`
      (sinon `make check` et tout run GPU abortent).
- [ ] Vérifier sur `src/snes/tutorials/ex19`:
  - `-dm_mat_type aijkokkos -dm_vec_type kokkos -pc_type gamg`
  - `-dm_mat_type aijcusparse -dm_vec_type cuda -pc_type hypre`
  - `-pc_type lu -pc_factor_mat_solver_type superlu_dist`
  - `-pc_type lu -pc_factor_mat_solver_type mumps`

## 4. Rodin avec PETSc + Kokkos

```sh
cmake -S rodin -B rodin/build-petsc \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=g++-14 \
  -DRODIN_USE_KOKKOS=ON \
  -DKokkos_DIR=$K/install-shared/lib/cmake/Kokkos \
  -DRODIN_USE_PETSC=ON \
  -DPETSc_DIR=$HOME/petsc-kokkos-install \
  -DRODIN_BUILD_UNIT_TESTS=ON \
  -DRODIN_BUILD_EXAMPLES=OFF \
  -DRODIN_INSTALL_RESOURCES=OFF \
  -DRODIN_SILENCE_EXCEPTIONS=ON
cmake --build rodin/build-petsc -j<N>
```

- `CMAKE_CXX_COMPILER` reste le compilateur hôte: `find_package(Kokkos)`
  installe `kokkos_launch_compiler`, seules les cibles liées à
  `Kokkos::kokkos` passent par `nvcc_wrapper` (voir `doc/Guides/Kokkos.dox`).
- `cmake/FindPETSc.cmake` lit `PETSc_DIR` ou la variable d'env `PETSC_DIR`.

- [ ] Configure + build.
- [ ] Revert de `SpatialComplexTest.cpp` si inutile avec gcc 14.
- [ ] `ctest --test-dir rodin/build-petsc/tests -R '^Kokkos_'` (29 tests).
- [ ] Tests PETSc et suite unitaire complète, noter les échecs restants.

## 5. Prototype Oseen GPU (hors repo)

Source: `~/Documents/oseen-gpu-proto` sur la machine locale (`proto.cpp`,
`CMakeLists.txt`, maillages `AtriumCoarse/Fine_FluidSolid_1mm.mesh`).
A copier sur demeter.

- Usage: `oseen mesh.mesh [steps] [rtol32] [rtol64] [64]`,
  `REFINE=n` pour n raffinements uniformes (x8 tets par niveau).
- Résultats locaux (RTX 5090, P1/P1, dt=1e-3, BiCGStab + Jacobi par blocs):
  - coarse fp32: 0.98 ms/step; fine fp32: 2.1 ms/step;
  - fine REFINE=1 (861k inconnues, 850 steps): fp32 32.1 ms/step
    (59 iters), fp64 37.9 ms/step (43 iters);
  - écart fp32 vs fp64 après 850 steps: L2 relative 0.67 en vitesse, 0.43 en
    pression (le flux de sortie reste a 3.4e-2).
- Le mur est le solveur, pas le GPU: les itérations croissent comme 1/h et
  avec le CFL, et fp32 demande plus d'itérations que fp64.

- [ ] Brancher le prototype (ou le module Kokkos de Rodin) sur PETSc:
      GAMG ou hypre BoomerAMG, ou fieldsplit vitesse/pression.
- [ ] Refaire l'étude de raffinement (REFINE=1, 2) avec multigrille.
- [ ] Comparer avec SuperLU_DIST GPU et MUMPS (CPU, MPI) en direct.
