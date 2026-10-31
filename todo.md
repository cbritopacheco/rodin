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
- [ ] Bug du prototype: BiCGStab sans détection de breakdown. Sur
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
- [ ] Brancher le prototype (ou le module Kokkos de Rodin) sur PETSc:
      GAMG ou hypre BoomerAMG, ou fieldsplit vitesse/pression.
- [ ] Refaire l'étude de raffinement (REFINE=1, 2) avec multigrille.
- [ ] Comparer avec SuperLU_DIST GPU et MUMPS (CPU, MPI) en direct.

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
- [ ] Mémoire locale: la matrice locale est dimensionnée a la compilation
      (8 sommets x 3 composantes par côté, soit 576 réels par thread pour
      un couple vectoriel). A revoir pour la performance (taille selon la
      géométrie, ou une équipe par cellule) au jalon benchmark.
- [x] `(Jac w).v` avec w vectoriel: fait au jalon 3.
- [ ] `Dot(A, Jac) : Dot(A, Jac)` (coefficient matriciel) reste a faire.
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
- [ ] Performance: dans `Mult<F, A>` le coefficient est évalué une fois par
      DOF local et par point de quadrature. Pour une GridFunction cela
      refait l'interpolation a chaque DOF; il faudra l'évaluer une seule
      fois par point de quadrature (cache dans le noyau) au jalon
      benchmark.
- [ ] `GridFunction` device construite par copie de `getData()` depuis
      l'hôte; l'interpolation ou la projection directement sur device
      reste a faire.

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
      sont rapatriées puis sommées dans le vecteur hôte, et
      `Default<Context::Kokkos<...>>` (un seul espace) choisit cet
      assemblage.
- [x] Intégrateurs négatifs (`lf -= ...`, `a - b`) déroulés par
      l'assemblage, linéaire et bilinéaire.
- [x] Tests `RodinKokkosVariationalLinearFormTest` comparés a l'hôte a
      1e-12 (scalaire, vectoriel, `VectorFunction`, `-=` avec `.over()`),
      et cas `a - b` ajouté a `RodinKokkosVariationalBilinearFormTest`.
- [ ] Performance: la somme des entrées dans le vecteur se fait sur
      l'hôte; une réduction atomique sur device (ou l'assemblage COO PETSc
      du jalon 6) évitera la copie.

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
      préassemblées et élimination des valeurs de Dirichlet sur l'hôte,
      identiques aux assemblages hôtes.
- [x] `Internal::forEachDeviceEntry` et `Internal::ConstrainedSystem`
      partagés par les assemblages de formes et de problèmes.
- [x] Tests `RodinKokkosVariationalProblemTest` comparés a l'hôte a 1e-12:
      Poisson avec lambda, constante, `.on()` et bord complet (triangles,
      mixte, tétraèdres), donnée `GridFunction`, problème vectoriel,
      Stokes mixte `(u, p, v, q)`, erreur sans incidence face-maille.
- [ ] Non supportés sur device: intégrateurs globaux, conditions
      périodiques et d'identification (`DirichletBC(u, A(v))`), qui
      lèvent une exception.
- [ ] `DirichletBC::on(1, 3)` avec des littéraux `int` ne compile pas
      sous nvcc (conversion rétrécissante dans `FlatSet<Attribute>{a1,
      a2, as...}`); les tests passent un `FlatSet<Geometry::Attribute>`.
- [ ] Performance: l'élimination et `setFromTriplets` restent sur l'hôte
      (jalon 6: COO PETSc).

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
- [ ] Formes PETSc préassemblées dans un problème device: lèvent une
      exception.
- [ ] MPI (`Context::Kokkos<MPI, B>`) non supporté.

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
