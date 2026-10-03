# Curved linear elasticity

This suite verifies isotropic linear elasticity on the fixed physical domain
$\Omega=\Phi((0,1)^d)$, for $d\in\{1,2,3\}$, where

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
-\operatorname{div}\sigma(u)=f,\qquad
\sigma(u)=\lambda\operatorname{div}(u)I+2\mu\varepsilon(u),\qquad
\varepsilon(u)=\tfrac12(Du+Du^T),\qquad u|_{\partial\Omega}=g.
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
 \|\varepsilon(u_h)-\varepsilon(u)\|_F^2\,\mathrm{d}x\right)^{1/2},
\qquad
E_\sigma=\left(\int_\Omega
 \|\sigma(u_h)-\sigma(u)\|_F^2\,\mathrm{d}x\right)^{1/2}.
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
$\|Au_h-b\|_2/\max(1,\|b\|_2)<10^{-11}$.

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
