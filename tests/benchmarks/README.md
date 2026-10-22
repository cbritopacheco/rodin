# Library performance benchmarks

This module measures the cost of library operations used by finite-element
assembly and solution. Its hierarchy follows the production modules, from
spatial arithmetic and reference data to mapped geometry, element kernels,
global assembly and solvers. Interpolation and projection are measured as
library operations in their own right.

The purpose is to locate work, distinguish its mathematical and implementation
dependencies, and evaluate proposed repairs without changing the discrete
problem. Local measurements are observations for a specified configuration,
not architecture-independent performance guarantees. This README therefore
defines measurement scopes and coverage gaps rather than recording timings,
speedup rankings or the chronology of individual investigations.

## Module architecture

```text
tests/benchmarks/
├── Workloads/                 Shared forms, inputs and numerical oracles
├── Rodin/
│   ├── Math/                  Spatial arithmetic and sparse algebra
│   ├── QF/                    Reference quadrature
│   ├── Geometry/              Connectivity, transformations and mapped points
│   ├── Variational/
│   │   ├── P1/                Low-order space-specific operations
│   │   └── H1/                High-order basis and local quadrature operations
│   ├── Solid/                 LinearElasticity element kernels
│   ├── Assembly/              Assembly stages and complete physical operators
│   ├── Solver/                Fixed-system solves and preconditioning
│   ├── Location/              Spatial search
│   ├── IO/                    Mesh input
│   ├── MPI/                   Distributed module specializations
│   └── PETSc/                 PETSc module specializations
└── convergence/               Verification-support profiling
```

Placement follows the operation being measured, not the physical context of
its caller. A spatial-matrix product belongs to Math; an element integral
belongs to Variational or Solid; global insertion belongs to Assembly or its
backend specialization. Shared workloads contain forms and reference checks,
not a second assembly policy. Qualified `Rodin/...` includes resolve to
production headers before mirrored benchmark fixtures.

Module and component executables are opt-in (`EXCLUDE_FROM_ALL`) and are not
registered with CTest. The compatibility executable `RodinBenchmarks` retains
the historical publication registration set. New standalone profilers are not
added to that publication workflow. Executables normally reside in
`build/tests/benchmarks/`; an explicitly configured runtime output directory
takes precedence.

The pressure-spectrum profiler under `convergence/` measures verification
support, not a production Stokes solver. A root-owned spectral calculation and
its result communication must not be interpreted as distributed assembly
imbalance. Such measurements do not count as production solver coverage.

## Sparse elimination

`RodinSparseEliminationBenchmarks` measures the public native sparse
`LinearSystem::eliminate` operation on P1 and degree-three stiffness-plus-mass
systems across Segment, Triangle, Quadrilateral, Tetrahedron, Hexahedron,
Pyramid and Wedge meshes. Sizes range from roughly one thousand to one hundred
thousand DOFs (the segment sweep also includes smaller systems). Random and
contiguous subsets prescribe 1%, 2%, 5%, 10%, 40%, 50%, 80%, or 90% of the
coefficients. The `placement` argument distinguishes random (`0`), contiguous
(`1`), and the complete assembled boundary (`2`, with `percent:0`). Boundary
indices come from `DirichletBC`; their values are coefficients of the exact
algebraic vector, not samples of an analytic boundary function.

The `Direct` family constructs directed upper-band CSC matrices of 1,000,
10,000 and 100,000 rows, with widths 8 and 64 and 1–10% scattered constraints.
The `missing:1` cases also omit columns and diagonals. They exercise the
absence of transpose storage and include diagonal insertion in the timed
operation. These are algebraic systems, not manufactured PDE convergence tests.

Assembly, subset construction, independent projector guards and matrix/load
restoration are excluded. FEM calls retain their compressed allocation. Direct
calls receive an exact-capacity copy so structural growth cannot disappear
through allocation reuse; copying is outside timing. All algorithm-owned
scratch and storage growth are inside timing. Manual timings report calling-
thread CPU seconds, including clock-read overhead; they are not wall-time or
whole-process throughput measurements. Each result is checked against the
matrix/load projector formula and the exact algebraic residual. Numerical guard
failures make the executable return failure.

The target is opt-in, serial, and not registered with CTest. For example:

```sh
cmake --build build --target RodinSparseEliminationBenchmarks
build/tests/benchmarks/Rodin/Math/RodinSparseEliminationBenchmarks \
  --benchmark_filter='SparseElimination/Direct/dofs:10000/width:64/percent:5/missing:1' \
  --benchmark_repetitions=7
```

Select a bounded family before running under a sampling profiler. Run one
benchmark process at a time. Use the same binary configuration, constraints,
coefficient restoration and timing scope for comparisons; retain repetition
and dispersion information on shared machines. The mathematical elimination
contract and algorithm are documented in the boundary-conditions guide and
in the sparse `LinearSystem::eliminate` Doxygen block.

## Mathematical measurement model

Let $\mathcal{T}_h$ be a mesh of $\Omega\subset\mathbb{R}^{D}$, with
reference cells $\widehat K\subset\mathbb{R}^{d}$ and maps
$F_K:\widehat K\to K$. The embedding dimension $D$ and reference dimension
$d$ are distinct. Let $V_h,W_h$ be finite-element trial and test spaces over
$\mathbb{F}\in\{\mathbb{R},\mathbb{C}\}$, with scalar, vector or matrix
values as supported by the operation. Their local basis counts are $N_u,N_v$.

For a quadrature rule $\{(\widehat x_q,w_q)\}_{q=1}^{N_q}$, a local
bilinear operator is represented by

$$
A^K_{ba}=\sum_{q=1}^{N_q}w_q\mu_K(\widehat x_q)
\mathcal{B}\bigl(\phi_a,\psi_b;c\bigr)(F_K(\widehat x_q)),
\qquad
\mu_K=\sqrt{\det(J_K^{\mathsf T}J_K)},\quad J_K=DF_K.
$$

Here $c$ denotes prescribed coefficients; $\mathcal{B}$ includes the required
pullbacks, derivatives and contractions. For full-dimensional cells the
measure factor is $|\det J_K|$. Point cells use the discrete point-measure
convention rather than this positive-dimensional derivative formula.

If $R_K^u,R_K^v$ restrict global coefficients to a cell, unconstrained assembly
has the algebraic structure

$$
A_h=\sum_{K\in\mathcal{T}_h}(R_K^v)^*A^K R_K^u.
$$

The adjoint $*$ is transpose-conjugation over $\mathbb{C}$. Constraint
expansion, insertion and completion are additional operations, not part of a
reference-basis microbenchmark. A component benchmark isolates one dependency
of this construction; an assembly consumer checks whether its cost matters to
the complete operation.

Comparisons record actual work, including $N_q,N_u,N_v$, component counts,
cells, global DOFs and stored nonzeros where applicable. A basis-value sweep
produces $N_qN_u$ scalar entries; a complete reference-derivative sweep produces
$dN_qN_u$. A general local matrix contains $N_vN_u$ entries. Equal mesh
parameters or polynomial degrees do not imply equal work across geometries.
Large dimensional differences require attribution to these counts and to
measured stages before they are explained as an architectural limitation.

## Physical contexts

| Context | Quantity measured and purpose | Important exclusions or gaps |
| --- | --- | --- |
| Mass, reaction and projection | Operator/load assembly and selected projection solves; distinguishes basis preparation, weighting and global insertion | Complete projection-solve/backend matrix, mixed/global couplings and component-storage attribution |
| Poisson, conductivity and reaction–diffusion | Diffusion/reaction operators and constrained manufactured problems; distinguishes geometry, coefficient and boundary-data evaluation | Vector/mixed constraints, broader boundary variants and distributed stage profiles |
| Complex Helmholtz | Complex assembly and fixed-system GMRES; separates operator construction from indefinite-system solution | Broader mapped/backend matrix and repeatable size-dependent solver attribution |
| Coupled reaction–diffusion | Block assembly and fixed-system solution; exposes coupling and scalar-promotion costs | Broader mixed spaces, backend parity and block-preconditioner stages |
| Stokes | Stable mixed velocity/pressure operators, pressure-mean coupling, local bindings and selected fixed-system solves | Wider degree/mesh/backend profiles; production solver timings remain distinct from the verification spectrum |
| LinearElasticity | Element and complete assembly kernels; measures symmetric-gradient/constitutive contractions | Higher degrees, complex and mapped/backend parity where supported |
| Nonlinear Poisson | Residual, tangent and source assembly at prescribed fields; isolates state-dependent kernels | Full Newton workflows, natural-boundary terms, PETSc/MPI parity and solver stages |

Darcy and NAFEMS remain deferred. Available mass benchmarks do not establish
completion of the remaining performance matrix.

### Coverage axes and applicability

The intended matrix distinguishes real/complex scalar fields; supported
scalar/vector/matrix ranges; P0, P0g, P1 and high-order H1 spaces; field degree;
reference-cell family; geometry degree and embedding; and native sequential,
OpenMP, PETSc-local and MPI execution. Reference-element drivers include H1
degree sweeps through six, while complete physical workloads generally cover
a smaller degree range. These are not interchangeable coverage claims.

The eight reference families are Point, Segment, Triangle, Quadrilateral,
Tetrahedron, Hexahedron, Pyramid and Wedge. Point algebra and basis checks do
not constitute a positive-dimensional PDE experiment. Constant spaces do not
provide conforming diffusion approximations; mixed problems require applicable
stable pairs. Unsupported combinations are excluded with a mathematical reason,
not silently counted as measurements or described as unimplemented physics.
A real-scalar PETSc installation does not certify complex PETSc behavior.

Mesh-based scaling requires at least three sizes; reference-element scaling
uses at least three quadrature orders. Field degree, geometry degree and
quadrature order remain independent inputs. Fixed quadrature defines a discrete
operator, not an exactness or convergence guarantee. In particular, rational
pyramidal bases and nonaffine maps require their own integration analysis.

## Measurement boundaries: what, how and why

Each benchmark defines its input state, timed operation and excluded work in
the owning fixture documentation. The following distinctions are essential.

| Scope | Measurement boundary | Interpretation |
| --- | --- | --- |
| Construction/lifecycle | Creation, copy, move, assignment or destruction of a named object | Separates setup/ownership costs from evaluation; construction is not implicitly included in a warm kernel |
| Reference preparation | Basis/rule construction or production at specified reference samples | Isolates degree, shape and quadrature dependence without global scatter |
| Mapped preparation | Coordinates, Jacobians, pullbacks, measure and cell binding | Exposes geometry-dependent work and distinguishes reference reuse from physical-state reuse |
| Local kernel | Preparation and/or entry evaluation for one cell, as explicitly named | Measures contractions before backend insertion; complete binding and repeated-entry reads are distinct scopes |
| Global assembly | Operator/load replacement, with setup/constraints/completion included only when stated | Determines whether an isolated change benefits its actual assembly consumer |
| Fixed-system solve | Solve on an already assembled operator, with preconditioner setup/reuse explicitly specified | Prevents assembly cost from being attributed to the solver or reuse from being confused with first setup |
| Interpolation/projection | Public DOF-functional interpolation or mass-system projection | Distinguishes coefficient construction from solving an approximation problem |

Cold/first-use, retained-reference and repeated-use states are separate
experiments. Changes of cell, rule, coefficient and geometry must remain visible
to the operation. Prepared-reference measurements do not represent ordinary
unprepared public calls. Setup, reference integration and numerical checks are
excluded from timing unless their inclusion is the stated measurement scope.

For distributed scopes, wall time is reported using the maximum participating
rank when the driver implements that protocol. Alignment and timing reduction
are excluded from the operation interval and documented separately. Owned and
halo work, empty shards and imbalance metadata are distinguished. No collective
is introduced into production code merely to simplify a local benchmark.
Strong scaling fixes global work; weak scaling fixes work per rank. Neither is
established by changing rank count and problem size without declaring which
quantity is held fixed.

## Supporting documents

[Measurement procedures](methodology.md) contain qualification, comparison and
execution instructions. The [scope inventory](coverage-ledger.md) tracks
available measurements and remaining work separately from the module design.
