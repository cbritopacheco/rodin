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

- [x] `g++-14 --version`, `gfortran-14 --version`: GCC 14.3.0
      (14.3.0-14ubuntu1). gcc 15 est le défaut sur demeter: toujours
      nommer `gcc-14`/`g++-14`/`gfortran-14` explicitement.
- [x] GPU et CUDA: RTX 2080 Ti, compute capability 7.5 (sm_75), driver
      595.99.02. CUDA 13.1 (V13.1.115) dans `/usr/local/cuda-13.1`, a mettre
      en tête du PATH. Ne pas utiliser `/usr/bin/nvcc` (CUDA 12.4, refuse
      gcc 14). `Kokkos_ARCH_TURING75=ON`, `--with-cuda-arch=75`.
- [x] MPI: OpenMPI 5.0.10 (apt, `/usr/bin/mpicc`, seul MPI; `/bin` est un
      lien vers `/usr/bin`). Pas CUDA-aware
      (`opal_built_with_cuda_support:value:false`).
- [x] OpenMPI 5.0.11 CUDA-aware dans `/opt/openmpi`, PETSc 3.26
      `arch-cuda-kokkos-ompi` (prefix `~/Documents/petsc-kokkos-ompi-install`),
      Boost 1.90 (`~/Documents/boost-1.90-ompi-install`), build Rodin
      `build-mpi` (`RODIN_USE_MPI=ON`).
- [x] CMake 4.2.3: exporter `CMAKE_POLICY_VERSION_MINIMUM=3.5`.
- [x] 112 coeurs, 750 Go de RAM: `-j64`.

## 2. Kokkos (build partagé, PIC)

PETSc en bibliothèque partagée exige un Kokkos partagé. Rodin doit utiliser
exactement le même Kokkos que PETSc (sinon deux runtimes Kokkos).

Kokkos 4.7.04 (tag `4.7.04`): CUDA 13 supporté depuis 4.7.01, fix de
compilation CUDA 13.1 dans 4.7.03 (kokkos#9020). Kokkos Kernels 4.7.04
correspondant pour PETSc.

```sh
K=$HOME/Documents/kokkos-src   # kokkos 4.7.04 cloné dans $K/kokkos
export PATH=/usr/local/cuda-13.1/bin:$PATH
export NVCC_WRAPPER_DEFAULT_COMPILER=g++-14
cmake -S $K/kokkos -B $K/build-shared \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=$K/kokkos/bin/nvcc_wrapper \
  -DCMAKE_CXX_STANDARD=20 \
  -DCMAKE_INSTALL_PREFIX=$K/install-shared \
  -DBUILD_SHARED_LIBS=ON \
  -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_SERIAL=ON \
  -DKokkos_ARCH_TURING75=ON \
  -DKokkos_ENABLE_IMPL_VIEW_LEGACY=ON \
  -DKokkos_ENABLE_COMPLEX_ALIGN=ON \
  -DKokkos_ENABLE_DEPRECATION_WARNINGS=ON
cmake --build $K/build-shared -j64 && cmake --install $K/build-shared
```

- [x] Build et install de Kokkos partagé: `$K/install-shared`
      (`libkokkoscore.so.4.7.4`), sans warning.
- [x] Optionnel: même chose en statique (`build`/`install`, sans
      `BUILD_SHARED_LIBS`, avec `CMAKE_POSITION_INDEPENDENT_CODE=ON`):
      `$K/install`.

## 3. PETSc 3.26 avec Kokkos, CUDA, hypre, SuperLU_DIST, MUMPS

PETSc 3.19 ne supporte pas Kokkos 4.7: prendre 3.26 avec Kokkos Kernels
4.7.04. PETSc n'a pas d'interface cuDSS (même en 3.26). PETSc 3.26.0
télécharge par défaut Kokkos et Kokkos Kernels 5.2.2 (minimum 4.3.00): le
`--download-kokkos-kernels-commit=4.7.04` est indispensable pour rester
aligné sur le Kokkos 4.7.04 partagé.

Source: `~/Documents/petsc-3.26` (tag `v3.26.0`).

```sh
export PATH=/usr/local/cuda-13.1/bin:$PATH
export OMPI_CC=gcc-14 OMPI_CXX=g++-14 OMPI_FC=gfortran-14
export CMAKE_POLICY_VERSION_MINIMUM=3.5
./configure PETSC_ARCH=arch-cuda-kokkos \
  --prefix=$HOME/Documents/petsc-kokkos-install \
  --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
  --with-cxx-dialect=20 --with-debugging=0 \
  COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3 \
  --with-cuda=1 --with-cuda-dir=/usr/local/cuda-13.1 --with-cuda-arch=75 \
  --with-kokkos-dir=$HOME/Documents/kokkos-src/install-shared \
  --download-kokkos-kernels --download-kokkos-kernels-commit=4.7.04 \
  --download-fblaslapack --download-scalapack --download-mumps \
  --download-metis --download-parmetis --download-superlu_dist \
  --download-hypre --with-make-np=64
make PETSC_DIR=$PWD PETSC_ARCH=arch-cuda-kokkos all
make PETSC_DIR=$PWD PETSC_ARCH=arch-cuda-kokkos install
```

- [x] Configure + build + install dans `~/Documents/petsc-kokkos-install`.
      Versions: Kokkos 4.7.4, Kokkos Kernels 4.7.4, hypre 3.2.0 (GPU, sans
      Umpire), SuperLU_DIST 9.2.1, MUMPS 5.9.1, METIS 5.1.0, ParMETIS 4.0.3,
      fblaslapack, ScaLAPACK. nvcc 13.1 avec `-ccbin mpicxx` (g++-14 via
      `OMPI_CXX`).
      Warning attendu: `mpi.mod` introuvable. Le `mpi.mod` d'OpenMPI apt est
      compilé par gfortran 15 et illisible par gfortran-14; PETSc se rabat
      sur `mpif.h`, sans conséquence.
- [x] MPI pas CUDA-aware: `export PETSC_OPTIONS="-use_gpu_aware_mpi 0"`
      pour tout run GPU. Attention, `make check` écrase `PETSC_OPTIONS`
      (makefile:133): passer l'option par
      `make ... EXTRA_OPTIONS="-use_gpu_aware_mpi 0" check`. Ainsi, tout le
      `make check` passe (ex19 1 et 2 procs, HYPRE/CUDA, CUDA, ex3k Kokkos
      Kernels, MUMPS, SuperLU_DIST, ex5f).
- [x] Vérifier sur `src/snes/tutorials/ex19` (2 procs, `-da_refine 3`):
  - `-dm_mat_type aijkokkos -dm_vec_type kokkos -pc_type gamg`: tourne sur
    GPU, mais KSP diverge (`DIVERGED_ITS`). Même divergence en `aij` CPU et
    en `aijcusparse`: c'est GAMG par défaut sur ce système non symétrique a
    4 composantes, pas Kokkos. Sans `-da_refine`, les trois types
    convergent identiquement (7 puis 6 itérations).
  - `-dm_mat_type aijcusparse -dm_vec_type cuda -pc_type hypre`: OK (3
    Newton, 18/24/28 itérations).
  - `-pc_type lu -pc_factor_mat_solver_type superlu_dist`: OK.
  - `-pc_type lu -pc_factor_mat_solver_type mumps`: OK.

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

- Prérequis apt installés: Boost 1.90, Eigen 3.4.0, HDF5 1.14.6, plus
  `git submodule update --init --recursive`.
- `ctest --test-dir build-petsc` ne trouve aucun test: utiliser
  `--test-dir build-petsc/tests`.

- [x] Configure + build (`-j36`): rc=0, 0 erreur, 224 warnings non bloquants
  (202 `-Wsign-compare`, 8 `-Wmaybe-uninitialized`, le reste divers).
- [x] Revert de `SpatialComplexTest.cpp`: inutile, le test compile et passe
  avec gcc 14.
- [x] `ctest --test-dir rodin/build-petsc/tests -R '^Kokkos_'`: 29/29.
- [x] Suite unitaire complète (dont 18 tests PETSc): 4958/4958, 1 skip
  attendu (`NamedFormConsistency.LinearElasticity/Segment`), 102 s en
  `-j36`. Aucun échec restant.

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

- Demeter (RTX 2080 Ti sm_75, CUDA 13.1, Kokkos 4.7.04 partagé, g++-14),
  build: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_CXX_COMPILER=g++-14 -DKokkos_DIR=$HOME/Documents/kokkos-src/install-shared/lib/cmake/Kokkos`
  (avec `NVCC_WRAPPER_DEFAULT_COMPILER=g++-14` et CUDA 13.1 dans le PATH).
  Temps de référence, `oseen <mesh> 100`, moyenne après warm-up, ms/step:

  | maillage (inconnues)   | fp32 step | iters | fp64 step | iters |
  |------------------------|-----------|-------|-----------|-------|
  | AtriumCoarse (41 416)  | 2.4-2.5   | 15    | 17.4-18.2 | 105   |
  | AtriumFine (114 448)   | 7.1-8.2   | 22-27 | 54.5-55.0 | 146   |

  Ms par itération: coarse 0.13 (fp32) et 0.16 (fp64), fine 0.22 et 0.35.
  Le fp64 est relativement cher sur Turing (débit fp64 au 1/32), et
  l'écart avec la 5090 est de x2.5 (coarse) a x3.4 (fine) en fp32.
- [x] Bug du prototype: BiCGStab sans détection de breakdown. Sur
      AtriumCoarse fp64 (rtol 1e-8), environ 1 run sur 3 tombe sur
      `(rh, r) = 0` a l'arrondi près (valeur exacte 1e-27 pour
      |rh||r| = 4e-10), `alpha = 0` puis `beta = rho/rhoOld` divise par 0
      a l'itération suivante: NaN dans x, propagé a tous les steps par le
      warm start (95 a 100 steps en échec). Non déterministe car l'ordre des
      `atomic_add` de l'assemblage change les arrondis. Pas un bug de
      réduction Kokkos (stress test de 600k réductions sans erreur). Fix
      proposé (non appliqué): redémarrer BiCGStab avec `rh = r`, `p = v = 0`
      quand `|rho| <= 1e-14 |rh| |r|` ou `omega == 0`. Validé sur 20 runs
      coarse (0 échec, même coût) et 3 runs fine.
      Clos: le prototype est hors du dépôt et remplacé par le module
      Kokkos de Rodin, dont le GMRES et le CG device (section 18) ne
      dépendent pas de ce BiCGStab.
- [x] Brancher le prototype (ou le module Kokkos de Rodin) sur PETSc:
      GAMG ou hypre BoomerAMG, ou fieldsplit vitesse/pression.
      Fait pour le module Kokkos de Rodin en section 18 (GAMG, hypre
      BoomerAMG via PETSc et en direct, fieldsplit sur Navier-Stokes).
- [x] Refaire l'étude de raffinement (REFINE=1, 2) avec multigrille.
      Clos: l'étude portait sur le prototype. Les mesures multigrilles
      du module Kokkos par taille de maillage sont en section 18.
- [x] Comparer avec SuperLU_DIST GPU et MUMPS (CPU, MPI) en direct.
      Fait en section 18, avec cuDSS.

## 6. DSL device (Kokkos/CUDA), jalon 1: formes bilinéaires P1 scalaires

Conception: artefact https://claude.ai/artifact/9F2T6hSwb8KHXrrxQEX7J1.
La syntaxe utilisateur est inchangée: même code `P1 fes(mesh)`,
`TrialFunction`, `Integral(Grad(u), Grad(v))`, `BilinearForm::assemble()`,
seul le type du maillage (`Mesh<Context::Kokkos<...>>`) change.

- [x] `Kokkos/Variational/Integral.h`: spécialisation
      `Integral<Dot<SFB<., DeviceP1>, SFB<., DeviceP1>>>` qui implémente
      `DeviceBilinearFormIntegrator::assemble()`: un kernel par intégrale,
      matrice locale dans un `Kokkos::Array`, triplets en mémoire device.
      Tables de quadrature aplaties par géométrie, filtrage par attributs.
      `integrate()`/`getPolytope()` lèvent une exception (pas d'évaluation
      hôte).
- [x] `Kokkos/Assembly/Kokkos.h`: `Assembly::Kokkos` sélectionné par
      `Assembly::Default<Context::Kokkos, Context::Kokkos>`, rapatrie les
      triplets et fait `setFromTriplets` (pas encore de COO PETSc).
- [x] `RealFunction` accepte un `KOKKOS_LAMBDA` sur `PointView`; son
      évaluation hôte lève une exception.
- [x] Ordre de quadrature aligné sur l'hôte, y compris le repli
      `k_tr + k_te - 2` des spécialisations P1 Grad.Grad et f Grad.Grad.
- [x] `Variational/P1/P1.h`: `static constexpr ElementType` devient
      `static const` (nvcc refuse le type non littéral). `P0.h` a le même
      motif, non touché.
- [x] Tests `RodinKokkosVariationalBilinearFormTest` (comparaison a l'hôte,
      tolérance relative 1e-12): masse, rigidité, somme, coefficient f,
      attributs, sur triangles, tétraèdres et maillage mixte quad/tri.
- [x] Règle d'ordre de l'hôte, reproduite sur le device: les handlers P1
      spécialisés avec dérivée (Grad.Grad, f Grad.Grad, Div(u).q, p.Div(v),
      Jac:Jac, f Jac:Jac, ...) prennent `k_u + k_v - 2`, ceux de type masse
      `k_u + k_v`, le handler générique la somme des ordres des feuilles. Le
      `-2` de Div(u).q n'est pas un bug: l'ordre 0 (point milieu) est exact
      sur simplexe pour P1 fois constante.

## 7. DSL device, jalon 2: P1 vectoriel et espaces mixtes

- [x] `Kokkos/Variational/P1.h`: `P1<SpatialVector<Real>, Mesh<Kokkos>>`,
      construit par `P1 fes(mesh, d)` comme sur l'hôte, même numérotation
      (composante r du sommet q en `q + r * nv`, local `a * d + r`).
- [x] `DeviceView`: valeur vectorielle, `Jacobian`, `Div`, somme, opposé et
      transposée de fonctions de forme, `f * A` et `A * f` (f scalaire),
      produit de Frobenius pour `Jac : Jac`.
- [x] `Integral` device généralisée a deux espaces P1 de rangs différents
      (essai vectoriel et test scalaire, et inversement): couplage
      vitesse/pression `Integral(Div(u), q)` et `Integral(p, Div(v))`.
- [x] Règle d'ordre `k_u + k_v - 2` étendue aux handlers hôte Div(u).q,
      p.Div(v), Jac:Jac, f Jac:Jac et (Jac f).v.
- [x] Tests `RodinKokkosVariationalVectorBilinearFormTest`: masse
      vectorielle, Jac:Jac, Div.Div, -Div.Div, gradient symétrique,
      coefficients, formes mixtes, sur triangles, tétraèdres déformés et
      maillage mixte quad/tri. Chaque test vérifie statiquement que
      l'intégrale device est bien sélectionnée.
- [x] Mémoire locale: la matrice locale est dimensionnée a la compilation
      (8 sommets x 3 composantes par côté, soit 576 réels par thread pour
      un couple vectoriel). Les matrices locales qui dépassent la capacité
      par thread sont intégrées avec une équipe par cellule (section 15).
- [x] `(Jac w).v` avec w vectoriel: fait au jalon 3.
- [x] `Dot(A, Jac) : Dot(A, Jac)` (coefficient matriciel): matrice
      constante et `Jacobian` d'une GridFunction comme coefficient, testés
      sur triangles, maillage mixte et tétraèdres. Les fonctions
      matricielles lambda device restent non supportées.
- [x] Jalons suivants: voir les sections 8 et suivantes.

## 8. DSL device, jalon 3: GridFunction et fonctions vectorielles

- [x] `DeviceView<GridFunction<P1 device, Data>>`: les coefficients sont
      copiés une fois dans une `Kokkos::View` a la construction de la vue
      (avant le noyau), puis interpolés par les fonctions de base P1 du
      point de quadrature. Scalaire et vectoriel, même numérotation que
      l'espace (`q + r * nv`).
- [x] `Grad(gf)` (P1 scalaire), `Jacobian(gf)` et `Div(gf)` (P1 vectoriel)
      sur device, via la jacobienne inverse du point.
- [x] `VectorFunction{f, g, ...}` sur device (au plus 3 composantes),
      composantes lambdas device, réels ou entiers; la CTAD accepte les
      lambdas `KOKKOS_LAMBDA` prenant un `PointView`.
- [x] Produits génériques `mult` (réel, vecteur, matrice) pour `f * A`,
      `A * f`, `M * v`, `M * N` entre fonctions et fonctions de forme.
- [x] `GridFunctionBase`: le chemin d'interpolation hôte n'est plus
      instancié sur un maillage Kokkos (trait `Context::IsKokkos`), ce qui
      permet d'écrire `GridFunction dg(dfes)` dans les expressions.
- [x] Tests `RodinKokkosVariationalCoefficientBilinearFormTest`: `g u.v`,
      `g Grad u.Grad v`, `(g + f) u.v`, `(g f) u.v`, convection
      `(Jac u) w.v`, Newton `(Jac w) u.v`, `Div(w) u.v`, `w Div(u).v`,
      `p Grad(g).v`, `VectorFunction` lambda et constante, comparés a
      l'hôte a 1e-12 sur triangles, maillage mixte et tétraèdres.
- [x] Performance de `Mult<F, A>`: le coefficient est lié une fois par
      point de quadrature (`bindBasis`, jalon 7).
- [x] Projection sur device des `GridFunction` d'un espace device
      (`gf = f`, `project(f)`, `project(region, f)`, avec un attribut ou un
      ensemble d'attributs): point de personnalisation
      `Variational::DeviceProjection<Context>` déclaré dans le coeur et
      spécialisé par le module Kokkos (`Kokkos/Variational/Projection.h`).
      La source (expression avec descente device, lambda `KOKKOS_LAMBDA`
      prenant un `PointView`, constante réelle ou vectorielle) est évaluée
      aux noeuds des DOFs des polytopes sélectionnés (cellules, faces, bord,
      interfaces); un DOF partagé prend la valeur du polytope d'indice
      maximal, comme le dernier visité sur l'hôte; les DOFs hors sélection
      gardent leur valeur. Tests `RodinKokkosVariationalProjectionTest`
      (P1, H1 degrés 2 et 3, scalaire et vectoriel) comparés a l'hôte a
      1e-12.
- [x] Projection device avec un prédicat: `project(region, f, pred)` avec
      un `KOKKOS_LAMBDA` sur `const Kokkos::Geometry::PolytopeView&`
      évalué dans le noyau de projection; un prédicat hôte sur
      `Geometry::Polytope` reste sur le chemin hôte.
- [x] `project(region, f, 2)`: les surcharges a prédicat sont contraintes,
      un attribut entier choisit la surcharge `Geometry::Attribute`.

## 9. DSL device, jalon 4: formes linéaires

- [x] `Integral` d'une fonction de forme test P1 device a valeur réelle:
      `Integral(v)`, `Integral(f, v)` (lowering en `Dot(f, v)`),
      `Integral(g, Div(v))`, `Integral(w, Grad(v))`,
      `Integral(Jacobian(w), Jacobian(v))`. Un thread par maille, vecteur
      local en registres, entrées `(indice, valeur)` écrites sur device.
- [x] `Internal::DeviceQuadrature`: tables de quadrature par géométrie de
      référence partagées par les intégrales bilinéaires et linéaires
      (même choix d'ordre que la `QuadratureRule` hôte).
- [x] `Assembly::Kokkos<Math::Vector<Real>, LinearForm<...>>`: les entrées
      sont sommées dans le vecteur sur le device (section 16), et
      `Default<Context::Kokkos<...>>` (un seul espace) choisit cet
      assemblage.
- [x] Intégrateurs négatifs (`lf -= ...`, `a - b`) déroulés par
      l'assemblage, linéaire et bilinéaire.
- [x] Tests `RodinKokkosVariationalLinearFormTest` comparés a l'hôte a
      1e-12 (scalaire, vectoriel, `VectorFunction`, `-=` avec `.over()`),
      et cas `a - b` ajouté a `RodinKokkosVariationalBilinearFormTest`.
- [x] Performance: la somme des entrées dans le vecteur se fait sur le
      device, dans l'ordre de l'hôte (section 16).

## 10. DSL device, jalon 5: problèmes et conditions de Dirichlet

- [x] `Assembly::Kokkos<IndexMap<Real>, DirichletBC<TrialFunction<P1 device>,
      FunctionBase<...>>>`: deux noyaux sur les faces de codimension un.
      Sélection par attribut (faces intérieures comprises, comme l'hôte)
      ou par incidence face-maille quand `.on()` est absent; chaque sommet
      contraint est attribué a sa première face (`atomic_min`), où la
      valeur est évaluée au sommet.
- [x] `Assembly::Kokkos<LinearSystem<SparseMatrix, Vector>, Problem<...>>`
      simple champ et multi-champs (`Problem(u, p, v, q)`, blocs par
      UUID): intégrateurs locaux et linéaires sur device, formes
      préassemblées, élimination des valeurs de Dirichlet et matrice
      compressée sur le device (`DeviceSystem`, section 16), identiques
      aux assemblages hôtes aux arrondis près (l'hôte OpenMP somme dans un
      ordre qui dépend de `OMP_NUM_THREADS`, écart de l'ordre de 1e-15).
- [x] `DeviceSystem` partagé par les assemblages de formes et de
      problèmes.
- [x] Tests `RodinKokkosVariationalProblemTest` comparés a l'hôte a 1e-12:
      Poisson avec lambda, constante, `.on()` et bord complet (triangles,
      mixte, tétraèdres), donnée `GridFunction`, problème vectoriel,
      Stokes mixte `(u, p, v, q)`, erreur sans incidence face-maille.
- [x] Intégrateurs globaux sur device: `Integral(Potential(k, u), v)`
      sur P1 scalaire ou vectoriel (noyau `Real(x, y)` ou
      `void(Matrix&, x, y)`), dans les formes, les problèmes simples et
      multi-champs et l'assemblage COO PETSc, attributs de test et de
      trial compris. Une paire de mailles par thread, dans l'ordre de
      l'hôte: Sauter-Schwab sur les triangles coïncidents, centroïdes
      sinon. Tests `RodinKokkosVariationalPotentialTest` et
      `Kokkos_PETSc_Problem.PotentialOnMixedMesh` comparés a l'hôte a
      1e-12 (triangles, tétraèdres, mixte, noyau singulier, négation).
- [x] Conditions périodiques et d'identification (`DirichletBC(u, A(v))`)
      sur device, Eigen (`DeviceSystem`) et PETSc (`DeviceCOOSystem`), via
      `DeviceConstraints`. Tests `RodinKokkosVariationalConstraintTest` et
      `Kokkos_PETSc_Problem.ConstraintsOnTriangles`.
- [x] Opérateurs `local - global`, `liste - global` et `corps - global`;
      `BilinearForm` accepte `=` pour les globaux, et `=` remplace
      désormais tous les intégrateurs (locaux et globaux).
- [x] `DirichletBC::on(1, 3)` et `over(1, 3)` avec des littéraux `int`:
      chaque attribut est converti explicitement en `Geometry::Attribute`.
- [x] Performance: l'élimination et `setFromTriplets` sont faits sur le
      device (section 16).

## 11. DSL device, jalon 6: assemblage COO PETSc

- [x] `Assembly::Kokkos<PETSc::Math::LinearSystem, Problem<...>>` simple
      champ et multi-champs (cible `Rodin::KokkosPETSc`, en-tête
      `Rodin/Kokkos/PETSc.h`): les triplets device des intégrateurs sont
      écrits dans une liste COO en mémoire device, avec l'élimination des
      valeurs de Dirichlet sur device (coefficients mis a zéro, relèvement
      par `atomic_add`, diagonale toujours listée, donc motif indépendant
      des conditions).
- [x] `MatSetPreallocationCOO` seulement si la matrice ou les coordonnées
      changent (comparaison sur device), puis `MatSetValuesCOO` avec un
      pointeur device; matrice `MATAIJKOKKOS` et vecteurs `VECKOKKOS` par
      défaut, second membre copié par `VecGetArrayWriteAndMemType`.
- [x] `PETSc::Variational::GridFunction` et `Problem` acceptent les
      contextes `Context::Kokkos<Local, B>` (trait `Context::IsLocal`),
      avec des données `VECKOKKOS`.
- [x] Tests `RodinKokkosPETScProblemTest`: système identique a
      l'assemblage device Eigen (Poisson, Stokes mixte), réassemblage sans
      nouvelle analyse du motif, cibles LHS et RHS, résolution CG Jacobi
      sur GPU comparée a SparseLU.
- [x] Formes PETSc préassemblées dans un problème device (simple champ et
      multi-champs): `BilinearForm` et `LinearForm` PETSc device autonomes
      (`MATSEQAIJKOKKOS`, `VECSEQKOKKOS`), entrées recopiées depuis les
      objets PETSc et combinées aux intégrateurs. Test
      `Kokkos_PETSc_Problem.PreassembledOnTriangles`.
- [x] MPI v1 (`Context::Kokkos<MPI, B>`, cible `Rodin::KokkosMPI`,
      `Rodin/Kokkos/MPI.h`): `Mesh<Kokkos<MPI, B>>` rend le shard d'un
      `Mesh<MPI>` résident (numérotation locale, fantômes compris) et marque
      les entités possédées dans la vue device; les intégrateurs ne visitent
      que les polytopes possédés. `P1<R, Mesh<Kokkos<MPI, B>>>`: membres hôte
      globaux (comme le P1 MPI hôte), membres device locaux. Problème PETSc
      simple champ assemblé sur le device en indices globaux
      (`MATMPIAIJKOKKOS`, COO `VecSetValuesCOO` pour le second membre),
      Dirichlet synchronisé par Vec fantôme (fixé partout dès qu'un
      processus le fixe). GridFunction PETSc en Vec Kokkos fantôme. Test
      `RodinKokkosMPIPoissonTest_np{1,2,4}` (build-mpi) contre l'hôte MPI.
- [ ] MPI: PBC, formes préassemblées, multi-champs, H1<2>, P0g et systèmes
      Eigen non distribués; vues device de `GridFunction<..., ::Vec>` non
      supportées.

## 12. DSL device, jalon 7: benchmark Oseen

- [x] Exemple `examples/Kokkos/OseenAssembly.cpp` (cible
      `KokkosOseenAssembly`, `KokkosOseenAssembly [mesh | n] [reps]`):
      même pas d'Oseen P1/P1 stabilisé assemblé dans PETSc par le DSL
      (`Problem` + COO device, un noyau par intégrateur) et par un noyau
      écrit a la main (matrice locale 16 x 16 par maille, intégrales P1
      exactes en forme close, sans quadrature). Résultats identiques a
      1e-15 près.
- [x] Optimisations génériques, temps DSL sur la grille 30 (162000
      tétraèdres, noyau manuel 6.4 ms dont 3.5 ms de `MatSetValuesCOO`):

      | étape | DSL | ratio |
      |---|---|---|
      | départ | | 38-58x |
      | cache COO, intégrateurs groupés par segment | 116 ms | 18x |
      | une jacobienne par point, coefficients sur le point | 94 ms | 14.5x |
      | fonctions liées une fois par point (`bindBasis`) | 84.7 ms | 13.1x |
      | bases P1 vectorielles creuses (`UnitVector`, `RowMatrix`) | 68.4 ms | 10.6x |
      | jacobienne affine réutilisée entre points (simplexes) | 59.0 ms | 9.1x |
      | bornes statiques (tétraèdres, dimension maximale) | 54.6 ms | 8.4x |
      | produits structurellement nuls sautés (`overlaps`) | 45.2 ms | 7.0x |

      Grille 16: 6.9 ms contre 1.1 ms (6.0x). Atrium `Coarse.mesh`
      (683704 tétraèdres): 234 ms contre 36.6 ms (6.4x).
- [x] Essai `TeamPolicy` (une équipe par maille): correct mais 2x plus
      lent (occupation faible a ~160 registres), abandonné.
- [ ] Objectif 10-20% non atteint, et pas atteignable avec une quadrature
      générique: le noyau manuel est une formule analytique, fusionne les
      8 termes en un seul passage et garde sa matrice locale en registres.
      Ecarts restants:
      - accumulateur local de 144 réels en mémoire locale (pile 2-8 Ko
        par thread, 160-250 registres), qui domine le coût par entrée;
      - 8 noyaux et 8 lectures-écritures du segment COO au lieu d'un;
      - quadrature d'ordre 2 (4 points) pour `Jacobian:Jacobian`, comme
        l'hôte, alors que l'intégrande est constant;
      - `IntegrationPointView` dimensionné pour 8 sommets même sur les
        tétraèdres.
- [ ] Pistes: fusionner les intégrateurs d'un même segment dans un seul
      noyau (somme des intégrandes), capacité de `IntegrationPointView`
      en paramètre, ordre exact pour les formes dérivées vectorielles
      (aussi côté hôte), découpage de la matrice locale par blocs de
      composantes.

## 13. DSL device, jalon 8: intégrales de bord, de faces et élasticité

- [x] `BoundaryIntegral`, `FaceIntegral` et `InterfaceIntegral` sur device
      (formes bilinéaires et linéaires), mêmes intégrandes que
      `Integral`: le noyau est générique en région
      (`Internal::DevicePolytopes`), parcourt les polytopes de dimension
      D - 1 et garde ceux de la région (bord: une maille incidente,
      interface: plus d'une) et des attributs demandés, comme l'hôte. La
      base P1 et les sommets sont ceux de la face, donc les mêmes degrés
      de liberté que l'hôte.
- [x] La région fait partie de `DeviceLayout`: les segments COO PETSc
      séparent les intégrales de mailles et de faces.
- [x] Sans incidence face-maille sur le maillage, une intégrale de bord ou
      d'interface lève une exception.
- [x] `DeviceView<BoundaryNormal>`: normale unitaire extérieure calculée
      sur device a partir de la jacobienne de la face, orientée par le
      barycentre de la maille incidente.
- [x] `LinearElasticityIntegrator` pour le P1 vectoriel device
      (`Kokkos/Variational/LinearElasticityIntegral.h`): somme des deux
      intégrales device `lambda Div(u) Div(v)` et
      `mu (Jac(u) + Jac(u)^T) : Jac(v)`, avec l'ordre de quadrature de
      l'hôte (`max(k_lambda, k_mu) + k_u + k_v`); le second terme
      s'accumule sur les entrées du premier.
- [x] Tests a 1e-12 contre l'hôte (`RodinKokkosVariationalFaceIntegralTest`):
      formes de bord, de faces et d'interface scalaires et vectorielles,
      avec attributs et coefficients, formes linéaires dont `Dot(n, v)`,
      élasticité avec Lamé constants et variables (triangles, mixte,
      tétraèdres), erreur sans faces. Problèmes Robin/Neumann et
      élasticité avec traction normale dans `RodinKokkosVariationalProblemTest`
      et `RodinKokkosPETScProblemTest`.
- [x] `.over()` des intégrales device renvoie le type device (surcharge
      dans chaque classe finale).
- [x] `Dot` de deux fonctions sur device (par exemple `Dot(Grad(g), n)`
      comme coefficient). `Grad`, `Jacobian` et `Div` d'une GridFunction
      évalués sur une face le sont dans la maille incidente, comme l'hôte:
      la seule sur le bord, celle du trace domain sur une interface (le
      noyau s'arrête s'il n'y en a pas). Testé sur le bord et sur les
      interfaces avec `traceOf`.
- [x] `Grad` et `Jacobian` d'une fonction test sur une face (bord et
      faces intérieures): gradient tangentiel de la face, comme l'hôte.
- [x] Élasticité sur les espaces H1 vectoriels device: une base commune
      `Kokkos::Variational::LinearElasticityIntegratorBase` et deux
      spécialisations (P1 et H1 vectoriels). Testé en H1 quadratique sur
      triangles et tétraèdres (`RodinKokkosVariationalH1Test`).

## 14. DSL device, jalon 9: exemples sur GPU

- [x] `examples/Kokkos/PDEs`: variantes device de `Poisson`,
      `MixedBoundary`, `Elasticity` (CG Eigen sur un `Problem` device) et
      `Seq_Poisson` (PETSc, hexaèdres 16^3, assemblage COO device). Le
      maillage hôte est construit puis copié sur le device; le reste du
      code est celui de l'exemple hôte.
- [x] Sorties comparées aux exemples hôte: maillages identiques,
      solutions à 1e-16 (Poisson, MixedBoundary, Seq_Poisson) et 3e-13
      (élasticité, CG à 1e-10).
- [x] `Mesh<Kokkos>::toHost()`: copie hôte du maillage (sommets,
      polytopes, attributs, transformations, incidences, nom), utilisée
      par `save()`.
- [x] `IO::XDMF` accepte un maillage et des `GridFunction` device: le
      maillage est recopié sur l'hôte par `write()` (à chaque instantané
      en `MeshPolicy::Transient`), le champ est recopié sur l'espace hôte
      de même numérotation (`P1::toHost(mesh)`). Tests
      `RodinKokkosIOXDMFTest` (fichiers identiques à l'hôte, maillage
      mobile) et `Kokkos_Mesh.ToHostRoundTrip`.
- [x] `DeviceView<Zero>` scalaire et vectoriel, pour `DirichletBC(u, Zero())`.
- [x] Correctif hôte: `Zero<SpatialVector>::getDimension()` manquait, et
      `VectorFunctionBase::getDimension()` se rappelait à l'infini.
- [x] Une directive `using namespace Rodin;` globale rend `Kokkos` ambigu
      dans les en-têtes Kokkos sous nvcc (`using namespace Kokkos;` de
      `Kokkos_CopyViews.hpp` résolu a l'instanciation). Inhérent au nom
      `Rodin::Kokkos`: documenté dans le guide Kokkos, la directive se met
      dans une fonction.
- [x] Portés: `BlockGravity` et `CantileverBeam`
      (`examples/Kokkos/Solid`, noyau device néo-hookéen
      `Kokkos/Variational/InternalVirtualWork.h`) et `NavierStokes`
      (`examples/Kokkos/PDEs/NavierStokes.cpp`). `Stokes` est porté au
      jalon 10.
- [x] `NavierStokes` device plus lent que l'hôte (1272 s contre 974 s):
      l'assemblage est 2 fois plus rapide sur le device, mais SparseLU
      domine. Pression passée de `P1` à `H1<1>` comme sur l'hôte: même
      espace, numérotation différente, et COLAMD remplissait davantage
      (60 s contre 44 s pour 100 pas).
- [x] Eigen n'est pas vectorisé dans une unité compilée par nvcc
      (`EIGEN_GPUCC` coupe `EIGEN_VECTORIZE`): SparseLU seul sur la même
      matrice, 9.7 s sous nvcc contre 7.9 s sous g++. `CG` et `SparseLU`
      creux sont maintenant instanciés dans `RodinSolver`, compilé par g++,
      avec des `extern template` activés par `RODIN_SOLVER_EXTERN_TEMPLATES`.
      La résolution sur le GPU passe par `Solver::CuDSS` (section 18).
- [x] `KSP::solve` et `SNES::solve` ignoraient les erreurs PETSc en
      Release (`assert` seulement): `-pc_factor_mat_solver_type cusparse`
      sur `seqaijkokkos` annonçait un succès en 1 ms avec une solution
      fausse. `*SetType`, `*SetFromOptions` et `*Solve` lèvent maintenant
      une `Alert::MemberFunctionException` avec le message PETSc.
- [x] `Darcy` reste vide côté hôte, rien a porter.
- [x] Lois de `Solid::InternalVirtualWork` sur device
      (`Kokkos/Variational/HyperElasticLaw.h`): Saint Venant-Kirchhoff,
      Mooney-Rivlin et Yeoh donnent la contrainte et le tenseur
      d'élasticité en forme fermée, `Kokkos::Variational::DeviceLaw` les
      associe aux lois hôte et des spécialisations génériques de la
      tangente (par thread ou par équipe) et du résidu contractent ce
      tenseur. Le néo-hookéen garde ses noyaux dédiés. Mêmes formes que
      l'hôte à 1e-12 et mêmes Newton pour les quatre lois. Sur hexaèdres
      cubiques, la tangente générique en équipe prend 0.29 à 0.32 s contre
      0.24 s pour le néo-hookéen. `HolzapfelOgden` et les lois actives
      lisent des entrées constitutives (fibres, activation), absentes du
      device: elles restent hôte.
- [x] Bogue cudafe++ (nvcc 13.1, gcc 14): `Integral(Grad(u), Grad(v)).over(1)`
      en premier terme d'une somme est mal parenthésé; le mettre ailleurs
      dans la somme. Contraintes nvcc gardées: noyaux en membres statiques
      publics, pas de première capture dans un `if constexpr` d'une lambda
      device, `using namespace` dans une fonction.

## 15. DSL device, jalon 10: H1 d'ordre 2, Taylor-Hood et P0g

- [x] `H1<2>` scalaire et vectoriel sur un maillage Kokkos simplicial
      (`Kokkos/Variational/H1.h`): degrés de liberté des sommets puis des
      arêtes, même numérotation que l'espace hôte construit sur
      `Mesh::toHost()` (`H1::toHost(mesh)`). Un maillage avec des
      quadrilatères lève une exception.
- [x] `P0g` sur device (`Kokkos/Variational/P0g.h`): un seul degré de
      liberté, base constante, pour le multiplicateur de moyenne de
      pression.
- [x] Degrés de liberté device génériques (`Kokkos/Variational/DOFs.h`):
      `VertexDOFs` (P1), `QuadraticDOFs` (H1 d'ordre 2) et `ConstantDOFs`
      (P0g, et espace d'essai implicite des formes linéaires). Les noyaux,
      les vues device et `DirichletBC` sont écrits contre le concept
      `DeviceSpace`; la capacité locale est fixée a la compilation par
      espace.
- [x] Ordre de quadrature identique à l'hôte, y compris pour les paires
      mixtes (H1 d'ordre 2 et P1, P0g).
- [x] `DirichletBC` sur H1 d'ordre 2: valeurs aux noeuds de Lagrange de
      la face (sommets et milieux d'arêtes, `getNode`).
- [x] Tests à 1e-12 contre l'hôte (`RodinKokkosVariationalH1Test`, 13
      tests): formes scalaires et vectorielles (masse, gradient,
      jacobienne, divergence, coefficients et `GridFunction` H1), couplages
      Taylor-Hood, formes de bord, de faces et d'interface, problèmes avec
      Dirichlet, Stokes Taylor-Hood avec multiplicateur P0g.
- [x] Correctif hôte: la forme linéaire `Integral` exige un intégrande
      scalaire (`ScalarValuedIntegrand`).
- [x] Un intégrateur multiplié par un réel (`viscosity * Integral(...)`)
      est accepté par l'assemblage device, Eigen et PETSc, comme la
      négation: le facteur multiplie les contributions.
- [x] `DeviceView` d'une `RealFunction` qui enveloppe une `RealFunction`,
      composante d'une `VectorFunction{RealFunction(f), ...}`.
- [x] `examples/Kokkos/PDEs/Stokes.cpp`: cavité entraînée Taylor-Hood
      P2/P1 avec multiplicateur P0g sur device, résolue par `SparseLU`.
      Sorties comparées à l'exemple hôte: maillage identique, vitesse à
      1e-14, pression à 2e-11 pour un maximum de 10.
- [x] H1 d'ordre K quelconque sur toutes les géométries sauf la pyramide
      (`LagrangeBasis`, `LagrangeDOFs`): base en monômes des coordonnées
      de référence, coefficients par inversion de la Vandermonde aux noeuds
      hôtes (Fekete, GLL) en précision étendue. La forme close de l'ordre 2
      reste sur les simplexes (Stokes identique au bit près). Numérotation
      à pas fixe par dimension, avec offsets sur les maillages mixtes.
      Tests à 1e-12: K=3 sur triangles et tétraèdres, K=2 sur quadrangles,
      hexaèdres et maillage mixte, K=3 sur quadrangles, refus des pyramides
      (29 tests).
- [x] Base de Lagrange des segments, quadrangles et hexaèdres évaluée en
      produit tensoriel des polynômes de Lagrange 1D aux noeuds GLL
      (poids barycentriques, sans branche), au lieu des monômes. Banc
      `KokkosH1_Assembly 30` (27000 hexaèdres Q2): Grad·Grad scalaire de
      880 à 575 ms, Jacobian·Jacobian vectoriel de 3100 à 2385 ms. Le
      reste du noyau par thread: géométrie 45 ms, bases 200 ms, produits
      300 ms (scalaire), la matrice locale débordant en mémoire locale.
- [x] Matrice locale au delà de 256 coefficients (et toujours au delà de
      `LocalCapacity`): une équipe par cellule sur device. Les points de
      quadrature sont traités par groupes en mémoire partagée (géométrie,
      valeurs des bases test et essai calculées en boucles séparées), puis
      chaque thread accumule des blocs de 8 coefficients (32 au delà de
      `LocalCapacity`). Le poids du polynôme 1D est sélectionné sans index
      par voie pour éviter la sérialisation de la mémoire constante. Banc
      `KokkosH1_Assembly`: scalaire Q2 (n=30) de 575 à 445 ms, scalaire Q3
      (n=14) de 3560 à 775 ms, vectoriel Q3 (n=14, 2197 hexaèdres, refusé
      avant) 5150 ms; le vectoriel Q2 reste par thread (2370 ms), ses
      valeurs matricielles ne tenant pas assez de points en mémoire
      partagée. Tests Q3 scalaire et vectoriel sur hexaèdres.
- [ ] Equipe vectorielle plus rapide: ne calculer par bloc que les lignes
      de bases utiles, éviter de recalculer géométrie et bases à chaque
      bloc de coefficients.
- [x] Tangente hyperélastique (`InternalVirtualWork`) au delà de
      `LocalCapacity`: une équipe par cellule, points et cinématique par
      groupes en mémoire partagée, gradients et tirés en arrière par
      paire point et fonction de base, blocs de 32 coefficients par thread.
      Test Q3 vectoriel sur hexaèdres (refusé avant) identique à l'hôte.
- [x] Gradient du déplacement de l'équipe hyperélastique sommé par paire
      point et fonction de base de l'état (addition atomique en mémoire
      partagée), au lieu d'un `interpolateJacobian` sériel par point.
      Tangente seule, déplacement nul: Q3 vectoriel n=14 (2197 hexaèdres)
      de 14,7 à 6,6 s; équipe forcée en Q2 n=20 de 7,2 à 3,2 s.
- [x] Equipe hyperélastique plus rapide que le thread dès P3 sur
      tétraèdres et Q2 sur hexaèdres. Une passe par groupe de points:
      chaque thread prend une fonction d'essai et un bloc de `TeamPairs`
      (4) fonctions test et accumule leurs blocs vdim x vdim en registres,
      en réutilisant les produits de la fonction d'essai. Gradients
      partagés entre espaces d'essai, test et état de même type, groupes
      de points deux fois plus grands (`TeamScratch` 47104 octets), noyau
      compilé pour 512 threads (`LaunchBounds`, 128 registres au lieu de
      254). RTX 2080 Ti, Neo-Hookean, Q2 vectoriel n=20: équipe de 3,2 s à
      0,46 s (thread 1,7 s); Q3 n=14: 6,6 s à 2,4 s. Seuil `LocalCapacity`
      abaissé de 8192 à 1024 entrées: P2 sur tétraèdres (900 entrées)
      reste par thread (n=30: 300 ms contre 455 ms en équipe, égalité sur
      RTX 5090), P3 sur tétraèdres (n=14) à égalité sur la 2080 Ti (180
      contre 190 ms). RTX 5090 (ares, sm_120), thread puis équipe: Q2 n=20
      1305 à 121 ms, Q2 n=30 1700 à 415 ms, P3 tétraèdres n=14 95 à 35 ms,
      Q3 n=14 en équipe 610 ms. Banc `KokkosH1_Assembly n répétitions
      degré [tetrahedra]` (degré 2 ou 3), tangente Neo-Hookean à petit
      déplacement.

## 16. Optimisations génériques de l'assemblage device

Mesures sur `KokkosSeq_Poisson 100` (hexaèdres Q1, un million de
cellules, CUDA sm_75).

- [x] Construction du maillage device: lecture directe de l'index des
      transformations du maillage local (une transformation absente est
      la transformation P1 par défaut, sans la construire), contrôles
      écrits sans vecteur intermédiaire par entité. Boucle sur les
      entités de 2.3 s à 0.15 s.
- [x] Copie de la connectivité: deux passes (tailles puis remplissage par
      index), remplissage direct du miroir hôte de la vue device
      (`upload` avec un remplisseur). De 2.83 s à 1.5 s; le reste est
      limité par les accès indirects de `getPolytope`.
- [x] COO PETSc: tri des coordonnées par ligne puis colonne sur le device
      (`sort_by_key`), coordonnées distinctes et valeurs sommées par
      coordonnée sur le device, de façon déterministe. PETSc reçoit des
      coordonnées triées et distinctes et ne trie plus. Premier
      assemblage de 6.96 s à 1.47 s (`MatSetPreallocationCOO` de 4.93 s
      à 0.62 s), copies hôte des coordonnées supprimées.
- [x] `IntegrationPointView::setPoint` réévalue le point en place au lieu
      de reconstruire la vue. Assemblages suivants de 525 ms à 423 ms
      (matrice de 416 à 334 ms, second membre de 92 à 72 ms).
- [x] Noyaux de points spécialisés par géométrie statique: quand le
      domaine d'intégration ne contient qu'une géométrie de la dimension
      de l'espace (triangle, quadrangle, tétraèdre ou hexaèdre), les
      noyaux évaluent la base de la transformation et la jacobienne à
      dimensions et nombre de nœuds statiques (`StaticBasis`,
      `StaticGeometry`). Les transformations identité et paramétriques
      stockées restent lues par entité. Noyau de la matrice de 337 à
      205 ms par appel, second membre de 73 à 4.3 ms. Hexaèdres: 125
      registres et 1.35 Ko de pile pour la matrice (matrice locale 8×8),
      94 registres sans pile pour le second membre.
- [ ] Ordre de quadrature de Grad·Grad sur hexaèdres Q1 (somme des degrés
      moins deux, comme l'hôte, soit 27 points): sur une transformation
      affine, le degré par direction est au plus 2 et la règle tensorielle
      d'ordre 2 (8 points) est exacte; sur une transformation trilinéaire
      l'intégrande est rationnel et aucune règle n'est exacte. Passer à un
      ordre par direction change le contrat numérique de l'hôte et du
      device (`doc/agents/quadrature.md`) et attend une décision. Gain
      attendu sur le noyau de la matrice: environ 3x.
- [x] Copie de la connectivité parallélisée côté hôte (`forEachOnHost`,
      OpenMP, lectures concurrentes du maillage local): copie de 1.79 s à
      environ 0.75 s avec 36 threads. Le reste est surtout le coût des
      fautes de page des tampons fraîchement alloués. Avec les 112 threads
      par défaut de demeter, la copie prend 1.1 à 1.4 s: l'attente active
      des threads ralentit les parties séquentielles.
- [x] Systèmes Eigen construits sur le device (`DeviceSystem`): les
      contributions restent là où les intégrateurs les écrivent, leur rang
      dans la liste de l'hôte est calculé à partir de leur position. Le
      tri par coefficient (puis par ligne pour l'élimination et le second
      membre) porte le numéro de la maille dans les bits de poids faible
      de la clé, ce qui laisse les contributions d'un coefficient dans
      l'ordre de l'hôte; un tri par insertion par groupe rattrape un tri
      instable. Les sommes suivent l'ordre de `setFromTriplets()`, les
      produits de l'élimination sont arrondis avant d'être sommés (pas de
      FMA), et seuls la matrice compressée et le second membre sont
      copiés sur l'hôte. Systèmes identiques bit à bit à l'assemblage
      précédent (Poisson, élasticité, bord mixte, Stokes). Poisson P1 sur
      706k tétraèdres: de 750 ms à 36.3 ms par assemblage; problème
      vectoriel: de 4516 ms à 92.6 ms. Il reste environ 10 ms de tris
      radix et 10 ms de copie vers l'hôte (mémoire non épinglée).
- [x] Assemblage COO PETSc (`DeviceCOOSystem`): l'élimination, le
      relèvement des valeurs affines des degrés de liberté identifiés et
      les intégrateurs linéaires listent leurs contributions au second
      membre à des positions fixes (comptage puis `parallel_scan`), au
      lieu de les ajouter par `atomic_add`. La liste est regroupée par
      ligne (`group` de `DeviceSystem`, le rang d'une contribution étant
      sa position) et chaque ligne est sommée dans l'ordre de la liste.
      Second membre identique bit à bit d'un assemblage à l'autre (test
      `Reassembly`). `KokkosSeq_Poisson 100`: noyaux du second membre de
      2.9 ms à environ 10 ms, premier assemblage inchangé (environ 1.1 s).
- [x] Lancements refusés signalés. Kokkos ne vérifie les lancements
      qu'avec ses contrôles de bornes: un noyau refusé au lancement (pas
      de mémoire pour ses tableaux locaux sur une carte presque pleine)
      laissait la matrice vide sans erreur. `KokkosH1_Assembly 40 1 2`
      sur la 2080 Ti: vectoriel Q2 en 89 ms avec 0 non nuls.
      `Kokkos::Internal::checkLaunches()` lève une exception après chaque
      `fill()` des assembleurs Kokkos et PETSc.

## 17. Temps de compilation des tests Kokkos

- [x] Tests Kokkos les plus lents découpés en plusieurs unités de
      compilation, liées dans le même exécutable. Dans build-mpi,
      `BilinearFormTest.cpp`, `VectorBilinearFormTest.cpp` et
      `CoefficientBilinearFormTest.cpp` (`tests/unit/Rodin/Kokkos/Variational`)
      passaient chacun plus de 26 minutes dans un seul nvcc en fin de build.
      Mesures dans build-petsc, compilés seuls en parallèle: 173 s, 635 s et
      398 s. Chaque unité paie environ 110 s d'en-têtes (un fichier d'un seul
      `TEST` sans noyau lourd: 112 s), et chaque fonction `check*` instancie
      ses noyaux; un fichier par `check` (12 unités) donne 320 s de temps mur
      mais double le temps CPU, et `BilinearFormTest.cpp` n'y gagne rien
      (172 s pour sa plus grosse fonction seule). Seuls les `check` lourds
      sortent dans leur propre fichier: `VectorBilinearFormFormsTest.cpp`,
      `VectorBilinearFormSymmetricGradientTest.cpp`,
      `VectorBilinearFormCoefficientTest.cpp` et
      `CoefficientBilinearFormVectorGridFunctionTest.cpp`. 7 unités, 349 s de
      temps mur au lieu de 635 s, la plus lente étant le gradient de
      `GridFunction` vectorielle (convection), qui ne se découpe plus.

## 18. Solveurs sur le GPU

Mesures sur demeter (RTX 2080 Ti, 32 coeurs), exemples
`examples/Kokkos/PDEs/Seq_*.cpp` pour PETSc.

- [x] `Solver::CuDSS` (`RODIN_USE_CUDSS`, cuDSS 1.0): LU, LDL^T ou LL^T sur
      le GPU, analyse et facteurs gardés comme `MUMPS`. La CSC de A est
      passée comme la CSR de A^T avec résolution transposée. Poisson 2D à
      1M inconnues: 6.6 s pour la première résolution contre 38 s pour
      `SparseLU`, puis 234 ms par refactorisation et 15 ms par résolution.
      Test `RodinCuDSSTest`.
- [x] `NavierStokes` Kokkos: problème et solveur gardés d'un pas à l'autre,
      `CuDSS` refactorise à chaque pas quand `RODIN_USE_CUDSS` est actif,
      `SparseLU` sinon. 160 pas (33 noeuds, Taylor-Hood avec P0g): 4.8 s
      au total, dont 1.7 s d'assemblage et 1.3 s de résolution, contre
      60 s avec `SparseLU` refait à chaque pas. Centre du couvercle
      identique à la référence sur les 160 pas.
- [x] Élasticité PETSc: CG avec GAMG et les modes rigides passés par
      `MatSetNearNullSpace` (construits à la main, la numérotation P1
      vectorielle est par composante): 28 itérations et 3.1 s à 822k
      inconnues, contre 14 s pour hypre, 27 s pour ILU, 19 s pour
      SuperLU_DIST; Jacobi ne converge pas.
- [x] Poisson 3D Q1 PETSc: Jacobi 208 ms à 1M inconnues, GAMG 6.1 s.
      hypre plantait à 1M (`CUDA ERROR out of memory`): BoomerAMG par défaut
      prend 5 Go à 262k. Avec `-pc_hypre_boomeramg_strong_threshold 0.7
      -pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_P_max 4`: 1.6 Go et
      647 ms à 262k (912 ms par défaut), 1.9 s à 1M. Documenté dans
      `Seq_Poisson.cpp`.
- [x] Navier-Stokes PETSc (`Seq_NavierStokes.cpp`, Taylor-Hood,
      `Problem::setFieldSplits`): FGMRES avec un split de Schur (`selfp`,
      GAMG sur la vitesse, quelques GMRES sur la pression), 6.1 s pour
      5 pas à 149k inconnues, contre 29 s pour ILU et 25 s pour
      SuperLU_DIST.
- [x] SuperLU_DIST sur un processus: l'offload GPU (actif par défaut) ne
      fait que 38 DGEMM, 12 ms sur 3 s de factorisation. Le mode
      `-mat_superlu_dist_3d -mat_superlu_dist_d 1` divise la factorisation
      par 3 (Poisson 3D, 33k: 3.1 s en 0.94 s; 262k: 130 s sur 32 coeurs
      en 2D, 48 s sur 8 coeurs en 3D, pendant un build).
      `SUPERLU_ACC_SOLVE=1` rend la factorisation 10 fois plus lente.
      PETSc n'a pas d'interface cuDSS.
- [x] `CuDSS` sur la CSR device sortie de l'assemblage Kokkos, sans
      passer par la matrice hôte: un `Problem` nommé avec
      `Math::LinearSystem<Kokkos::SparseMatrix<Real, Backend>,
      Math::Vector<Real>>` garde sa matrice sur le device (CSR, clés
      ligne-majeure, mêmes sommes que la CSC hôte donc mêmes coefficients au
      bit près, tableaux réutilisés d'un assemblage à l'autre), seul le
      second membre est copié. `CuDSS` lit les tableaux en place (mode
      NORMAL, `cudssMatrixSetCsrPointers` à chaque refactorisation), via le
      trait `Internal::CuDSSOperator`. Pas de forme bilinéaire préassemblée.
      Test `RodinKokkosVariationalDeviceSystemTest` (Poisson, périodique,
      Stokes mixte, CuDSS). `NavierStokes` 160 pas: 3.3 s au lieu de
      3.6 s (assemblage 1.19 s en 1.06 s, résolution 1.06 s en 0.94 s).
- [x] Krylov sur Kokkos-Kernels sans PETSc (`RODIN_USE_KOKKOSKERNELS`,
      Kokkos-Kernels 4.7.4 sur le même Kokkos, ordinal et offset `int`):
      `Solver::CG` et `Solver::GMRES` spécialisés pour un système à
      `Kokkos::SparseMatrix`, Jacobi, CSR lue en place. CG écrit sur Kokkos
      (sémantique du CG Eigen, mise à jour fusionnée avec deux réductions,
      `SPMVHandle` gardé pendant la résolution: 208 ms sans, 121 ms avec),
      GMRES de Kokkos-Kernels avec redémarrage. Poisson 3D Q1 à 1M
      (`Krylov_Poisson.cpp`): 171 ms et 158 itérations à 1e-10 contre
      210 ms et 159 itérations pour PETSc Jacobi, assemblage 211 ms contre
      1056 ms. Test `RodinKokkosSolverKrylovTest` (CG et GMRES contre
      `SparseLU`).
- [x] hypre appelé directement (`RODIN_USE_HYPRE`, `FindHYPRE.cmake`,
      hypre 3.2.0 CUDA de petsc-kokkos-install, AmgX non installé):
      `Solver::BoomerAMG`, PCG de hypre préconditionné par un V-cycle
      BoomerAMG, pour un système à `Kokkos::SparseMatrix`. CSR passée par
      `HYPRE_IJMatrixSetValues2` en mémoire device, réglages GPU de PETSc
      (PMIS, ext+i, l1-Jacobi, élimination de Gauss au niveau grossier),
      seuil, niveaux agressifs et P_max réglables, norme 2 non
      préconditionnée. MPI et hypre initialisés au besoin, finalisation
      laissée à l'application comme `MUMPS`. hypre fait toujours au moins
      une itération et ne signale
      pas la convergence pour b = 0 (traité avant). Poisson 3D Q1
      (`BoomerAMG_Poisson.cpp`, 0.7, 2, 4): 644 ms et 42 itérations à 1M
      contre 1.9 s et 32 itérations pour PETSc en norme non préconditionnée
      (27 en norme préconditionnée), 404 ms contre 592 ms à 262k. L'écart
      d'itérations n'est pas expliqué. Test `RodinKokkosSolverBoomerAMGTest`.
- [x] Hiérarchie BoomerAMG gardée d'une résolution à l'autre, avec la
      copie hypre de la matrice dont elle vient, comme
      `KSPSetReusePreconditioner`: le CG tourne sur une copie de la matrice
      courante et le V-cycle sur la matrice gardée
      (`HYPRE_PCGSetPrecondMatrix`, setup du préconditionneur vide), donc
      une hiérarchie périmée ne coûte que des itérations. `setup()` la
      refait, `clear()` la libère, les réglages de la hiérarchie la
      libèrent, une copie du solveur ne la partage pas, un changement de
      taille lève une erreur. Poisson 3D Q1 (0.7, 2, 4): setup 463 ms puis
      260 à 280 ms par résolution à 1M, contre 644 ms quand tout était
      refait; 343 ms (initialisation MPI et hypre comprise) puis 90 à
      113 ms à 262k. Test: matrice doublée résolue avec l'ancienne
      hiérarchie, `setup`, `clear`, copie.
