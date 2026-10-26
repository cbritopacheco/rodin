# Atrial wall contraction demo

This example simulates transient, finite-strain contraction of a material-4
atrial wedge wall. A prescribed activation plane drives a fiber-reinforced
Neo-Hookean material. Reference-domain inertia is integrated with Newmark
average acceleration, using Rodin's `Integral`, `InternalVirtualWork`,
`SparseLU`, and `NewtonSolver`.
Activation is prescribed; no electrophysiology equation or evolving
Hill–Maxwell internal state is solved.

## A plane wave starting at the anatomical base

In the original `Coarse.mesh` volume, label 6 identifies the basal cap. The
original tetrahedra have volume material 1. In the wedge-only wall, the cap
is absent and its boundary label is preserved on the adjoining rim faces.
The wave direction is determined by the basal rim.

The inner vertices of rim 6 define a least-squares plane through their mean
position $b$. Its unit normal $d$ is oriented toward the wall centroid.
The projected coordinate and activation are

$$
s(X)=d\cdot(X-b),\qquad
 a(X,t)=\exp\!\left[-\left(\frac{s(X)-vt}{w}\right)^2\right],
\qquad c(X,t)=-c_{\max}a(X,t).
$$

The wave crest starts on the fitted basal plane at $t=0$ and moves in the
positive $d$-direction. The initial displacement and velocity are zero;
the initial acceleration balances the basal active force. The last crest
position is three widths past the maximum wall projection.

This pulse propagates only in the positive $d$-direction. Points with
$s(X)<0$ receive only the initial Gaussian tail. Complete propagation along
the wall requires a different activation construction, such as a wall-distance
wave.

The fitted origin, inward normal, projected coordinate range, and landmark
information are written to `<stem>.geometry.json`. `WaveCoordinate` is a
static nodal XDMF field for checking the plane orientation.

## Simple landmark-based fibers

The default left-appendage apex landmark is

$$
X_{\mathrm{LAA}}=(46.39638367,187.44307667,-205.03317767)\ \mathrm{mm}.
$$

Coordinates are used instead of vertex or triangle numbers because MMG
renumbers entities. An apex patch is selected using inner-triangle centroids
within 3 mm of this point, together with the corresponding outer triangles.
A finite patch avoids imposing a point-only Dirichlet condition in the wall
volume. Its radius and center can be changed from the command line.

A single harmonic coordinate is solved with Rodin:

$$
-\Delta\psi=0\quad\text{in }\Omega_{\mathrm{wall}},\qquad
\psi=0\quad\text{on rim 6},\qquad
\psi=1\quad\text{on the appendage patch},
$$

with zero normal flux on the remaining boundary. A temporary boundary label
8 is used for this solve and then restored, preserving all original wall and
rim labels.

For each wedge, the known extrusion vectors supply a unit transmural
normal $m$. The harmonic gradient is projected into the wall tangent
plane and the fiber is constructed as

$$
\ell=\frac{(I-mm^T)\nabla\psi}{\|(I-mm^T)\nabla\psi\|},
\qquad f=\ell\times m.
$$

Thus the fibers are tangent to the wall and follow contours of the
base-to-appendage coordinate. They are stored as one constant unit direction
per wedge. There is no transmural rotation in this one-layer model.

This is a simplified use of the Laplace-gradient frame in Piersanti et al.,
*Modeling cardiac muscle fibers in ventricular and atrial electrophysiology
simulations*, [arXiv:2101.10960v1](https://arxiv.org/abs/2101.10960v1), equations
(3) and (5). It does not implement their regional atrial bundle selection,
pulmonary-vein rules, or a validated anatomical fiber reconstruction.

Critical points are handled explicitly. If the harmonic tangent vanishes,
the projected direction toward the apex is used; at a remaining pole, a
Cartesian tangent is chosen. `FiberFallback` flags these cells: 0 means the
harmonic construction, 1 the geometric-apex fallback, and 2 the pole fallback.
`FiberCoordinate`, `LongitudinalDirection`, `WallNormal`, and
`FiberDirection` are exported for inspection.

## Transient mechanics

On the reference wall $\Omega_0$, the displacement $u$ satisfies

$$
\rho_0\ddot u-\mathrm{Div}_X P(F,a)=0,
\qquad F=I+\nabla_X u.
$$

The weak residual is

$$
R(u;w)=\int_{\Omega_0}\rho_0\ddot u\cdot w\,dX
       +\int_{\Omega_0}P(F,a):\nabla_Xw\,dX.
$$

The stress combines a compressible Neo-Hookean material with the static
constitutive branch of `ActiveFiberLaw`. The prescribed active target is
$c=-c_{\max}a$. “Static” here describes the constitutive branch: displacement,
velocity, and acceleration evolve through the momentum equation. Active
internal-variable kinetics, chamber pressure, and fluid coupling are absent.
All rim faces with labels 1, 2, 3, 5, and 6 have fixed displacement; the inner
and outer wall faces, labeled 4 and 7, are traction-free.

### Time integration and initial conditions

Newmark average acceleration uses $\beta_N=1/4$ and $\gamma_N=1/2$. With
$\Delta t$ the time step, the predictors are

$$
u^*=u_n+\Delta t\,v_n+\Delta t^2(1/2-\beta_N)\ddot u_n,
\qquad
v^*=v_n+\Delta t(1-\gamma_N)\ddot u_n.
$$

Newton's method solves for $u_{n+1}$ with

$$
\ddot u_{n+1}=\frac{u_{n+1}-u^*}{\beta_N\Delta t^2},
\qquad
v_{n+1}=v^*+\gamma_N\Delta t\,\ddot u_{n+1}.
$$

The tangent contains the elastic tangent and the consistent mass matrix
scaled by $1/(\beta_N\Delta t^2)$. The last accepted displacement is the
initial Newton iterate; the kinematic predictor enters the inertial residual.
No physical damping or algorithmic damping is added.

The initial state has $u_0=v_0=0$. Because the activation pulse is already
present at the base, its consistent acceleration is computed from

$$
M\ddot u_0=-f_{\mathrm{int}}(u_0,a(0))
$$

on the free degrees of freedom. Displacement, velocity, and acceleration
vanish at every fixed rim. The first XDMF state records these initial data;
each subsequent state follows one nonlinear time step.

The pulse traverses the wall over
$T=(\max_X s(X)+3w)/v$. `--steps N` sets $\Delta t=T/N$ and writes $N+1$
states including the initial data. The default $N=24$ is a coarse smoke test;
it does not establish time accuracy for inertial oscillations. Time-step
refinement is required for quantitative dynamics. Stability properties of the
linear Newmark scheme do not guarantee nonlinear convergence or temporal
resolution for this wall.

### Parameters and units

Coordinates are in mm, time in seconds, and stresses in kPa. `--density`
accepts reference mass density in $\mathrm{kg/m^3}$. Its internal value is

$$
\rho_{\mathrm{internal}}=10^{-9}\rho_{\mathrm{SI}},
\qquad
1060\,\mathrm{kg/m^3}
=1.06\times10^{-6}\,\mathrm{kPa\,s^2/mm^2}.
$$

| Parameter | Default |
|---|---:|
| Passive Young's modulus | $50\,\mathrm{kPa}$ |
| Passive Poisson ratio | $0.3$ |
| Active fiber stiffness | $20\,\mathrm{kPa}$ |
| Peak target Green–Lagrange strain magnitude | $0.03$ |
| Reference mass density | $1060\,\mathrm{kg/m^3}$ |
| Wave speed | $200\,\mathrm{mm/s}$ |
| Gaussian width | $6\,\mathrm{mm}$ |
| Basal boundary label | 6 |
| Appendage patch radius | $3\,\mathrm{mm}$ |
| Time intervals | 24 |

Mechanical quadrature order 4 is an example policy for the non-polynomial
constitutive law. Quantitative use requires spatial, quadrature, and temporal
refinement and material calibration.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target AtriumContraction -j4
build/examples/Heart/AtriumContraction/AtriumContraction \
  --mesh resources/examples/Heart/Atrium/AtriumWall_Repaired_Wedges_1mm.mesh \
  --output output/atrium-transient/Contraction \
  --base-label 6 --apex 46.39638367 187.44307667 -205.03317767 \
  --apex-radius 3 --steps 24 --contraction 0.06 --density 1060 \
  --speed 200 --width 6
```

`--young`, `--poisson`, and `--active-stiffness` configure the material.
`--stop-after` limits the number of positive-time steps while retaining the
initial state and the time step set by `--steps`. `--fibers-only` exports the
geometry and fiber fields without solving mechanics; its zero displacement,
velocity, and acceleration are inspection placeholders. `--help` lists the
options.

The preflight check finds the exact minimum of each six-node wedge's
Jacobian determinant over its reference cell. The determinant is affine in
the triangle coordinates and quadratic through the thickness, so evaluating
the triangle vertices and thickness-quadratic minima is sufficient. A
nonpositive reference determinant rejects the input. The same test is applied
to each converged deformed cell before saving a mechanics snapshot.

## Geometry and organization

[`AtriumContraction.cpp`](AtriumContraction.cpp) contains the geometry checks,
basal-plane fit, fiber-coordinate problem, constitutive setup, and Newmark
time loop. Its CMake target links Rodin's Solid, Variational, IO, and Solver
modules. The mesh resources are in `resources/examples/Heart/Atrium/`:

- `Coarse.mesh`: source atrial volume and labeled boundary.
- `AtriumWall_Repaired_Surface.mesh`: wall surface prepared with a local
  handle-to-disk repair and MMGS adaptation.
- `AtriumWall_Repaired_Wedges_1mm.mesh`: one outward layer of six-node wedges,
  with a nodal extrusion distance of $1\,\mathrm{mm}$.

The wall construction uses the material-4 surface. Inlet and outlet regions
are open, with their labels transferred to rim quadrilaterals. The local
surface repair changes the source topology; the resulting geometry is a
modeling assumption requiring anatomical assessment.

## ParaView

Open `<stem>.xdmf`, keeping its companion HDF5 files in the same directory.
Apply **Warp By Vector**, choose `Displacement`, and use scale factor 1 for
physical displacement in mm. Color by `Activation` and play the time series.

For fibers, apply **Cell Centers**, then line glyphs oriented by the cell
vector `FiberDirection`. Color the wall by `FiberCoordinate` and inspect
`FiberFallback` near the apex patch. `ActiveExtension` stores the prescribed
target strain. `Velocity` and `Acceleration` are nodal fields in
$\mathrm{mm/s}$ and $\mathrm{mm/s^2}$. The CSV records their maximum magnitudes,
time, maximum displacement, activation maximum, Newton iterations, final
residual, and minimum scaled Jacobian. Geometry metadata also records density,
time step, and Newmark parameters.
