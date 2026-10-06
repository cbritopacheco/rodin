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
