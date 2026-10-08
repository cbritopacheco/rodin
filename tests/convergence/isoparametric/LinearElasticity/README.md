# Curved linear elasticity

This suite verifies isotropic linear elasticity on the fixed physical domain
$\Omega=\Phi((0,1)^d)$, for $d\in\lbrace 1,2,3\rbrace$, where

$$
\Phi(\xi)=\xi+0.1\xi_0^2e_{d-1}.
$$

The quadratic geometry is represented exactly on every refinement level,
including boundary entities and MPI halos. Its Jacobian determinant is one
in two and three dimensions; in one dimension it is $1+0.2\xi_0>0$.
Thus domain approximation does not change during this field-refinement study.
The map checks in the existing curved-geometry suite provide an independent
geometry control.

## Formulation and manufactured data

With fixed dimensionless Lamé coefficients $\lambda=1.5$ and $\mu=0.5$,
the equations and full essential boundary condition are

$$
-\mathrm{div}\sigma(u)=f,\qquad
\sigma(u)=\lambda\mathrm{div}(u)I+2\mu\varepsilon(u),\qquad
\varepsilon(u)=\tfrac12(Du+Du^T),\qquad u\rvert_{\partial\Omega}=g.
$$

The smooth field is $u_i(x)=(i+1)\exp(s)$, with
$s=\sum_{j=0}^{d-1}x_j$. The independently specified source is

$$
f_i(x)=-\exp(s)\left[\mu d(i+1)
 +(\lambda+\mu)\sum_{j=0}^{d-1}(j+1)\right].
$$

Both volumetric and shear terms are present. This is a compressible,
displacement-only problem, not a nearly incompressible or locking test.
The analytic tensor references are specified independently of the discrete
Jacobian-to-tensor calculation:

$$
\varepsilon_{ij}(u)=\tfrac12(i+j+2)\exp(s),\qquad
\sigma_{ij}(u)=\left[\lambda\tfrac{d(d+1)}2\delta_{ij}
 +\mu(i+j+2)\right]\exp(s).
$$

## Measurements and acceptance

In addition to displacement L2 and H1-seminorm errors, the suite integrates

$$
E_\varepsilon=\left(\int_\Omega
 \Vert\varepsilon(u_h)-\varepsilon(u)\Vert_F^2\thinspace \mathrm{d}x\right)^{1/2},
\qquad
E_\sigma=\left(\int_\Omega
 \Vert\sigma(u_h)-\sigma(u)\Vert_F^2\thinspace \mathrm{d}x\right)^{1/2}.
$$

The expected orders are $p+1$ for displacement L2 and $p$ for its
H1 seminorm, strain L2 and stress L2, assuming the requisite solution and
dual regularity. Every adjacent interval checks finite positive errors,
strict decrease, and an order within $0.55$ of $p+1$ for displacement L2
or $0.45$ of $p$ for derivative quantities. These finite-resolution windows
are acceptance policies, not rigorous error bounds.

| Field degree | Grid points per axis | Geometry degree | Interpretation |
| --- | --- | --- | --- |
| P1 | $5\to9\to17$ | P2 | Superparametric field refinement |
| P2 | $3\to5\to9$ | P2 | Isoparametric field refinement |

Here $h=1/(n-1)$. All seven positive-dimensional geometries are exercised:
segment, triangle, quadrilateral, tetrahedron, pyramid, hexahedron and wedge.
There is no spatial elasticity convergence problem in zero dimensions.
Native Eigen and real-PETSc local/MPI builds share the same acceptance logic;
MPI registrations cover one through four ranks. MPI norm integration counts
owned physical cells only before a global squared-norm reduction.
A separate constant-vector norm check compares against
$\sqrt{d|\Omega|}$, where $|\Omega|=1.1$ in one dimension and one otherwise.

Separate constant P1 and physical-affine P2 patches check errors below
$10^{-9}$. A physical affine field pulls back to degree two under $\Phi$;
it is not asserted to belong to the P1 mapped space. A physical quadratic
field need not belong to P2 under this map and is not used as an exact patch.

Assembly uses quadrature order 11 and norm integration uses order 13.
Separate sensitivity solves increase assembly/norm orders to 16/18 or
tighten the CG relative tolerance from $10^{-13}$ to $10^{-14}$; every error
quantity must change by less than $10^{-6}$ relatively. The independently
computed coefficient residual satisfies
$\Vert Au_h-b\Vert_2/\max(1,\Vert b\Vert_2)<10^{-11}$.

An omitted-volumetric-term solve retains the original manufactured forcing
and trace. It must produce at least twice the baseline displacement L2,
strain L2 and stress L2 errors. This negative control demonstrates that
these observables reject an incorrect operator. It does not assert the
incorrect solution has a particular convergence order.

The tests are labelled `convergence;slow` (plus `petsc` and `distributed`
where applicable). Each geometry is independently registered with a
30-minute safety timeout for instrumented CI; this is not a performance claim.
Pyramid registrations share a resource lock to avoid overlapping the largest
curved workloads within a CTest run.

## Exact-domain vector lift and metric controls

The shared `LiftedErrorNorm` also accepts real and complex vector fields.
Let $x=\Phi(\xi)$ and $x_h=\Phi_h(\xi)$, and define
$u_h^\ell(x)=u_h(x_h)$. With component-row Jacobians, the chain rule is

$$
J(u_h^\ell)=Ju_h(x_h)D\Phi_h(\xi)D\Phi(\xi)^{-1}.
$$

The derivative defects of the field, geometry and total errors are evaluated
in one owned-cell quadrature traversal. An observer applies the linear strain
and stress laws to these defects, before reducing squared tensor norms.
Pulling back an already formed stress tensor is not the same operation.
Reference and represented cells are paired by logical indices and ordered
vertices, not by coordinate tolerances.

An independent metric oracle takes the represented map to be the identity
and the exact map to be $\Phi(\xi)=\xi+a\sin(\pi\xi_0)e_{d-1}$,
with $a=0.1$. The physical affine field is $u(x)=\mathbf{1}+Ax$, where

$$
A_{ij}=(i+1)(j+1)+\delta_{i0}\delta_{j,d-1},\qquad
c=Ae_{d-1},\qquad C^2=\Vert c\Vert_2^2,\qquad s=a\pi.
$$

The field error vanishes for a P2 solve on the identity mesh. Geometry and
total errors therefore coincide. Their displacement norms are

$$
E_{G,0}=aC/\sqrt{2},\qquad
E_{G,1}=\begin{cases}
C\sqrt{(1-s^2)^{-1/2}-1},&d=1,\cr
sC/\sqrt{2},&d\in\lbrace 2,3\rbrace.
\end{cases}
$$

For $d\in\lbrace 2,3\rbrace$, their constitutive norms are independently given by

$$
E_{G,\varepsilon}=s\sqrt{(C^2+c_0^2)/4},\qquad
E_{G,\sigma}=\frac{s}{\sqrt{2}}
\sqrt{2\mu^2 C^2+(d\lambda^2+4\lambda\mu+2\mu^2)c_0^2}.
$$

In one dimension these reduce to $E_{G,\varepsilon}=E_{G,1}$ and
$E_{G,\sigma}=(\lambda+2\mu)E_{G,1}$. The nonsymmetric matrix has
distinct relevant row and column norms in higher dimensions; the oracle
therefore detects an incorrectly transposed or left-multiplied lift.
These tests use $n=3$, quadrature order 18 and absolute tolerance $10^{-9}$,
on all seven geometries, in native and real-PETSc local/MPI configurations.

A separate complex-vector interpolation oracle uses $(1+\mathrm{i})u$.
Its displacement norms must be $\sqrt{2}$ times the real values, with
vanishing field error. It checks complex vector interpolation and lifting
on local and distributed meshes; it does not certify a complex-PETSc
elasticity solve.

## Smooth fields on approximated geometry

The sine map above is interpolated with degree-two geometry on each mesh,
defining $\Omega_h=\Phi_h((0,1)^d)$. The exponential manufactured field,
forcing and trace are evaluated in physical coordinates on $\Omega_h$;
the continuous represented-domain solution is therefore known independently
of the geometry approximation. This is not a boundary-data approximation test.

For $x=\Phi(\xi)$ and $x_h=\Phi_h(\xi)$, define the vector defects

$$
e_F(x)=u_h(x_h)-u(x_h),\qquad
e_G(x)=u(x_h)-u(x),\qquad e_T=e_F+e_G.
$$

All lifted norms use the exact-domain measure
$\det D\Phi(\xi)\thinspace \mathrm{d}\xi$. The strain and stress errors apply
their linear constitutive laws to $Je_F$, $Je_G$ and $Je_T$;
they are not obtained by composing a represented-domain tensor with the lift.
Both triangle and reverse-triangle inequalities are checked for all four
quantities with an absolute roundoff allowance of $10^{-11}$.

Subject to regular maps, the stated approximation and dual-regularity
assumptions, degree-two geometry predicts the following algebraic orders:

| Component | Displacement L2 | Displacement H1 seminorm, strain L2 and stress L2 |
| --- | --- | --- |
| Represented-domain and lifted field errors | $p+1$ | $p$ |
| Geometry defect | $3$ | $2$ |
| Total exact-domain error | $\min(p,2)+1$ | $\min(p,2)$ |

The P1 hierarchy is $n=5,9,17$; P2 uses $n=3,5,9$, except the
segment uses $n=5,9,17,33$ to exclude its coarse pre-asymptotic regime.
All adjacent intervals require finite positive, decreasing errors, with the
same $0.55$ displacement-L2 and $0.45$ derivative-rate margins used above.
These observations certify the specified finite hierarchies, not a uniform
asymptotic theorem for arbitrary meshes.

At $n=5$, separate controls increase only the assembly order $11\to16$,
only the norm order $13\to18$, or tighten only the solver tolerance
$10^{-13}\to10^{-14}$. Every represented and lifted error quantity must
change by less than $10^{-6}$ relatively. An omitted-volumetric-term P2
solve retains the correct manufactured forcing and trace. Its represented,
lifted-field and total errors must exceed twice their correct counterparts;
the geometry defect must remain identical, since it does not depend on the
discrete operator. All seven geometries and the native, real-PETSc local/MPI
and OpenMP configurations use the same observables and acceptance policies.
The approximated-domain cases have separate `Approximated` CTest registrations,
retaining the slow-test timeout and pyramid resource lock.

## Linear and cubic geometry with a representable displacement

The sine-map studies additionally use geometry degree
$q\in\lbrace 1,3\rbrace$ and field degree $p=\max(2,q)$.
The asymmetric affine field $u_\ast(x)=\mathbf{1}+Ax$, with $A$
defined above, has constant strain and stress and zero body force.
Its represented-domain pullback belongs to the field family. Consequently,
represented and lifted field errors can be required to remain below
$10^{-9}$, independently of the nonzero geometry defect.

For $x=\Phi(\xi)$ and $x_h=\Phi_h(\xi)$, define

$$
e_G(x)=A\bigl(\Phi_h(\xi)-\Phi(\xi)\bigr),\qquad
J_G=A\bigl(D\Phi_h(\xi)D\Phi(\xi)^{-1}-I\bigr).
$$

The geometry strain and stress defects are

$$
\varepsilon_G=\tfrac12(J_G+J_G^T),\qquad
\sigma_G=\lambda\mathrm{tr}(J_G)I+2\mu\varepsilon_G.
$$

Their exact-domain norms use the same measure and component-row chain
rule as the quadratic studies. The prescribed matrix has a nonzero last
column and first entry, so the displacement and constitutive observables
are sensitive to the chosen sine-map perturbation; arbitrary matrices
could cancel a particular observable. The total and geometry norms must
agree within the $10^{-9}$ field-reproduction budget, separately for
displacement L2, displacement H1 seminorm, strain L2 and stress L2.

| Geometry degree $q$ | Field degree $p$ | Grid points $n$ | Nominal displacement L2 / derivative orders |
| --- | --- | --- | --- |
| 1 | 2 | $3\to5\to9$; segment $5\to9\to17$ | $2/1$ |
| 3 | 3 | $3\to5\to9$; segment $5\to9\to17$ | $4/3$ |

Each adjacent interval requires finite, positive, decreasing geometry
and total errors. Displacement L2 rates must lie within $0.55$ of
$q+1$; displacement-gradient, strain and stress rates must lie within
$0.45$ of $q$. These checks measure the nonvanishing interpolation defects
of this finite experiment; an approximation upper bound alone does not
prove these two-sided rate windows.

At $n=5$, assembly order $11\to16$, norm order $13\to18$ and solver
tolerance $10^{-13}\to10^{-14}$ are varied independently. Both
represented and lifted field reproduction remain required in all four
solves. Every geometry and total observable must change relatively by
less than $10^{-6}$. The independent asymmetric metric and
omitted-volumetric-term controls above remain separate.

The segment regularity argument for these linear/cubic sine interpolants
is given in the [complex Helmholtz geometry methodology](../Helmholtz/README.md#linear-and-cubic-geometry-with-a-representable-complex-field).
In one dimension the represented physical interval is unchanged; the
nonzero lifted geometry defects measure the prescribed parametrization
comparison. Full essential traces and positive Lamé coefficients retain
the displacement-only coercive formulation on regular represented domains.

All seven positive-dimensional geometries have degree-specific entries
for native and real-PETSc local/MPI ranks one through four. Sequential
and OpenMP configurations are distinct verification gates. These entries
retain slow labels, 1800-second watchdogs, MPI processor counts and the
shared pyramid resource lock. The registration specifies the intended
matrix, not evidence of a completed execution.

## Mixed traction on exact quadratic geometry

The real-PETSc boundary target uses the same exact quadratic domain and
Lamé coefficients, with $\Gamma_D=\Phi(\lbrace \xi_0=0\rbrace)$ and
$\Gamma_N=\partial\Omega\setminus\Gamma_D$. Reference-face attributes
are assigned before partitioning and mapping. The manufactured physical
traction is $t=\sigma(u_\ast)n$, with the mapped outward unit normal $n$.
The weak problem is

$$
\int_\Omega\bigl(\lambda\mathrm{div}u_h\mathrm{div}v_h
+2\mu\varepsilon(u_h):\varepsilon(v_h)\bigr)\thinspace dx
=\int_\Omega f\cdot v_h\thinspace dx+\int_{\Gamma_N}t\cdot v_h\thinspace ds,
\qquad v_h\rvert_{\Gamma_D}=0,
$$

with $u_h=I_h^{\Gamma_D}u_\ast$ on the essential boundary. The nonempty
essential face removes rigid motions; Korn's inequality and $\mu>0$ give
coercivity. This is a compressible mixed-boundary study, not a uniform
locking-free assertion. The existing flat nearly-incompressible studies
remain separate and unchanged.

The exponential physical field above has P1/P2/P3 studies. P1 uses
$n=5,9,17$, except Tetrahedron retains the finer $n=9,17,33$ path from
the flat mixed-traction fixture. P2/P3 use $n=3,5,9$. Both adjacent
intervals must decrease, with displacement L2/H1-seminorm rate floors
$1.65/0.75$, $2.45/1.55$, and $3.45/2.55$, respectively. The same
regularity and finite-resolution qualifications apply as above.

At $n=3$, physical-affine and physical-quadratic patches use field degrees
two and four. Their pullbacks under the quadratic map have degrees at most
two and four, respectively. The asymmetric affine field is
$u_\ast=\mathbf{1}+Ax$, with $A$ defined in the vector-lift section above;
the quadratic field is $u_{*,i}=1+(i+1)(\sum_jx_j)^2$.
Both displacement errors must be below $10^{-9}$.

Two separate affine-P2 controls retain the original manufactured data.
One omits only the traction load; the other omits only the volumetric
bilinear term, retaining the full physical traction. The correct patch
must meet the absolute reproduction budget, whereas each wrong operator
must give $E_0>10^{-3}$ and $E_1>10^{-2}$. These are dimensionless
separation budgets on the prescribed unit-scale domain, not index-matching
tolerances. Entity correspondence and boundary attributes remain logical.

Assembly order $16\to18$, norm order $18\to20$, and CG relative
tolerance $10^{-13}\to10^{-14}$ are varied independently on the smooth
P2 field; both errors must change relatively by less than $10^{-6}$.
CG/Jacobi requires a positive convergence reason and an independently
recomputed coefficient residual below $10^{-11}$.

The target registers all seven geometries locally and at MPI ranks one
through four, separately under sequential/OpenMP assembly. Real-PETSc
storage and the existing vector workload are reused; no solver or library
implementation is added. Geometry installation and physical data evaluation
remain local; global assembly, solution and owned-cell norm sums require
communicator participation. Registrations retain slow labels and the shared
pyramid lock. Tetrahedral traction groups reserve a two-hour execution budget
for the resolved P1 hierarchy; the other geometries retain 30-minute budgets.
CI separates local and each MPI rank count for tetrahedral traction, so the
cost of five complete groups is not accumulated in one job. These are
scheduling policies, not numerical error budgets or performance bounds.
A registered test is not, by itself, evidence of completed numerical validation.
