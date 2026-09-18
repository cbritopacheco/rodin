# Kelvin ball

This directory separates evaluation of the Kelvin-ball objective from the
shape update.

`KelvinBallSphere` is the validation executable. It constructs and discretizes
the spherical one-twenty-fourth chamber, computes the resistance coefficients
with a stabilized P1--P1 Stokes discretization and rotational Nitsche coupling,
then reconstructs the complete fluid mesh from the 24 proper cubic rotations
and solves the full-domain problems independently. The chamber and full-domain
values of

\[
  \rho = \frac{|c|}{\sqrt{kq}}
\]

are printed together with their relative difference. The reconstructed mesh
is written to `KelvinBallSphere.xdmf`.

```sh
KelvinBallSphere --n=13 --outer-radius=2
KelvinBallSphere --h=0.125 --outer-radius=2 --save-mesh
```

`KelvinBallSphere --help` lists the shared geometry and discretization options.
The dimensionless Nitsche penalty and pressure-stabilization factor may be set
with `--penalty` and `--stabilization`; their defaults are 320 and 0.05.

The chamber calculation solves the three translational and three rotational
states directly on the fundamental chamber. The cut meshes need not match:
quadrature points are rotated to the paired face and located by an
attribute-filtered AABB tree. No pressure gauge is added to this coupled
family, since the rotational identification removes its constant modes. The
independent full-domain states each use one global mean-pressure constraint.

The mesh attributes are:

- `2`: body;
- `3`: fluid;
- `7`: outer container;
- `13`: body--fluid interface;
- `31`, `32`, `33`, and `35`: chamber cuts.

`KelvinBallOptimization` contains the optimization loop. It calls the same chamber
metric evaluator used by `KelvinBallSphere`, constructs the Hadamard density, projects
the descent direction onto the null space of the volume differential, and
normalizes the complete direction in the nodal infinity norm. There is no line
search.

The source is organized by mathematical responsibility. `Sphere` owns the
background chamber and its MMG discretization, `Metrics` owns the stabilized
Stokes families and resistance contraction, `RotatedNitscheIntegrator` owns
the AABB trace map and all mapped-face assembly kernels, and `SewedOutput` owns
the 24-copy reconstruction and field transport. `KelvinBallOptimization` owns
the staged optimization procedure. Namespace-level operations in `Common` are
limited to data and mesh preparation shared by these classes.

The objective and volume differentials are first identified with their Hilbert
gradients. With the \(\rho\)-gradient `rhoGradient` and the volume gradient
`volumeGradient`, the code evaluates

\[
  \lambda^*=-\frac{DV[\nabla\rho]}{DV[\nabla V]},
  \qquad
  \xi_\rho=\nabla\rho+\lambda^*\nabla V.
\]

The coefficient is evaluated through the assembled volume differential rather
than a separately assembled approximation of the Hilbert metric. Consequently,
the discrete direction satisfies \(DV[\xi_\rho]=0\) before normalization.

When Rodin is configured with `RODIN_MULTITHREADED=ON`, variational forms are
assembled with OpenMP. When `RODIN_USE_UMFPACK=ON`, the chamber and gradient
identification systems use UMFPACK; otherwise they use Eigen SparseLU. The
selected backends are printed at startup and recorded in `kelvin-ball.csv`.

The optimiser constructs the one-twenty-fourth chamber from a tetrahedral
`UniformGrid`. The initial level set is

\[
  \phi_0(x)=\lvert x\rvert-1,
\]

so its negative region is the spherical part of the chamber. The chamber outer
radius is \(L=2\) by default and may be set with `--outer-radius=L`; sewing the
24 chambers produces the cube \([-L,L]^3\). The resolution is specified either
by the number of grid points on an edge or by a requested mesh size. In the
latter case, the resolution is rounded up and the effective mesh size is
`L / (n - 1)`.

```sh
KelvinBall --n=17 --iterations=20
KelvinBall --h=0.125 --iterations=20
KelvinBall --outer-radius=3 --h=0.1666666667 --iterations=2
```

The shape density is extended with the regularized metric

\[
  \int_D \ell^2 \nabla\theta:\nabla w + \theta\cdot w\,dx,
  \qquad \ell=4h.
\]

The dimensionless multiplier can be changed with, for example,
`--regularization=6`; this does not alter the Stokes or MMG parameters.

Running `KelvinBall --help` prints every option and its default value. With no
arguments, the executable uses `--n=13`, `--outer-radius=2`, `--iterations=1`,
`--penalty=320`, `--stabilization=0.05`, `--regularization=4`, and `--step=0.1`.

Each iterate is written to two XDMF series. `KelvinBall.xdmf` contains the
chamber design and chamber states. On iterations followed by an update, its
`Advected` field is the transported distance on the exact mesh passed to MMG.
`KelvinBallSewed.xdmf` contains the complete
24-copy design with the sewn distance and deformation fields, together with a
complete fluid grid carrying the six sewn velocities and six sewn pressures.
It also contains the fully sewn `Advected` field. `KelvinBallMMG.xdmf` records
the reconstructed and subsequently optimized MMG mesh, before the next
finite-element spaces are constructed.
The console and `kelvin-ball.csv` also report the three side lengths of the
axis-aligned bounding boxes of the chamber mesh and the fully sewn mesh.
The CSV records one row per evaluated design, including the run parameters,
mesh and material counts, minimum, mean, and maximum tetrahedron mean-ratio
quality, mean element size, MMG reconstruction parameters, resistance metrics,
volume error, update prediction and realization, normalized direction
derivatives, extension residuals and jumps, and the wall-clock duration of each
numbered stage. Fields
that do not exist for the initial design or a `--state-only` run are `nan`.

A one-iterate output smoke test is

```sh
KelvinBall --n=9 --iterations=1 --step=0.01
```

Two iterations are required to exercise one geometry update and recompute the
state on the resulting design:

```sh
KelvinBall --n=13 --iterations=2 --regularization=4 --step=0.01
```

The transported-distance quadrature defaults to order 8 and can be changed
with `--advection-quadrature=<order>`. Order 2 exactly integrates only the
zero-displacement P1 mass-product limit. The composition with a nonzero
characteristic displacement is not polynomial on one source element and
requires a higher-order rule in practice.

The uniform grid supplies the background resolution. MMG discretizes the
initial spherical level set and each subsequently advected level set. The
outer boundary and all chamber cuts are marked as required geometry; their
labels and planarity are checked after every reconstruction. Since the body
interface crosses the chamber cuts, their triangulations may nevertheless
change locally. The rotational Nitsche terms locate paired traces by AABB
search and do not require matching cut meshes.

During level-set transport, a characteristic crossing an identified chamber
cut is continued through its rotated partner. The scalar projection weakly
matches the two nonconforming traces. These operations express the topology of
the sewn chamber; the outer boundary remains a physical boundary. The Eikonal
distance is first projected into the same weakly matched trace space, so the
subsequent advection increment does not include an unrelated trace correction.

The MMG parameters are tied to the effective background size: `hmin=0.1 h`,
`hmax=10 h`, `hausd=0.1 h^2`, and gradation `2`. Automatic angle detection
is disabled because the fixed chamber faces are protected explicitly. The
element count before and after each reconstruction is reported. Each
reconstruction is followed by
an MMG mesh-optimization pass with the same size and geometry parameters. The
fixed chamber faces remain required during this pass. The reconstructed update
is not identical to the advected deformation, so its realised objective and
volume changes are reported against their shape-derivative predictions at the
following iterate.
