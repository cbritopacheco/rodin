# Numerical contracts

This page collects cross-cutting numerical rules that otherwise hide in
module docs. Read it before changing interpolation, projection, assembly,
solvers, quadrature, or residual/tangent code.

## Interpolation, projection, and coefficients

`GridFunction::operator=(FunctionBase)` applies the finite element space's
degree-of-freedom functionals to the expression. This is interpolation in the
space's DOF sense, not an `L2`-orthogonal projection.

For an `L2` projection, assemble and solve the mass problem:

```cpp
TrialFunction u(Vh);
TestFunction  v(Vh);
Problem l2(u, v);
l2 = Integral(u, v) - Integral(f, v);
Solver::CG(l2).solve();
```

Do not read coefficients as values except where the space guarantees that
contract. P1 nodal spaces store vertex values. High-order H1 uses
Fekete/Dubiner machinery, so coefficients are not plain nodal values.

### Distributed interpolation

PETSc region interpolation on an MPI mesh uses the backend-independent
MPI specialization of `Interpolation<FES>`. The space supplies its context,
scalar field, DOF mappings, and functionals. For each owned global DOF $i$,
an eligible incident entity $K_i$ is selected by its smallest distributed
entity index.
The existing finite-element pullback and functional define the coefficient:

$$
c_i = \ell_{i,K_i}(f).
$$

The predicate and source function must agree on replicas of the same entity.
The vertex overlap supplies the incident entities for owned DOFs in $P_0$,
$P_1$, and $H^1$. For globally supported $P_{0g}$, the smallest eligible entity
is selected globally and its evaluated coefficients are communicated to the
DOF owner, which may have an empty shard. For a nonconstant source in
$P_{0g}$ this defines a deterministic interpolation sample, not a mean or an
$L^2$ projection. Interpolation of a representable constant is independent of
the selected source.

All selected coefficients are evaluated before the destination is updated.
Only owners commit coefficients, followed by an owner-to-ghost refresh.
Unselected coefficients are preserved, including when the eligible region is
empty. Direct mutable grid-function access retains its separate `acquire()` /
`flush()` contract; it is not an indexed remote-write interface. $P_0$
interpolation requires cells; its DOFs do not define face interpolation.

## Problem sign convention

Rodin problem assignment states a residual equation with everything moved to
the left:

```cpp
problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
```

The assembler routes bilinear terms to the operator and linear terms to the
right-hand side with the appropriate sign. For nonlinear Newton forms, build
the tangent and residual so the assembled system is

```text
J(x) dx = -F(x)
```

Do not negate the residual twice.

## Real and complex form convention

Rodin writes forms in trial-first, test-second order. `Math::dot(lhs, rhs)`
conjugates `rhs`, so a complex-valued form has the discrete actions

```text
L(v)   = v* b
a(u,v) = v* A u
```

Here `*` denotes the conjugate transpose. A `LinearForm` is therefore
conjugate-linear in its test argument over the complex numbers. A
`BilinearForm` is sesquilinear: linear in the trial argument and
conjugate-linear in the test argument. Over the reals these reduce to the
ordinary linear and bilinear cases; the class names cover both scalar domains.

Assembly stores `b[i] = L(psi_i)` and `A[i,j] = a(phi_j, psi_i)`. The local
integrator returns this entry directly. Do not add a backend-level conjugation:
`Integral(c * u, v)` represents `c * u * conj(v)`, not
`conj(c) * conj(u) * v`. Eigen evaluates the action as `v.dot(A * u)` and
PETSc as `VecDot(A * u, v)`. Eigen conjugates its first operand; PETSc
conjugates its second operand. Their argument order therefore differs even
though both evaluate the same action. Likewise, PETSc linear-form action
uses `VecDot(b, v)`. Do not infer one backend's dot-product convention from
the other backend's API spelling.

Sesquilinearity is not Hermitian symmetry. Whether `A == A*` is a separate
property of the particular weak form and coefficients.

## LinearSystem lifetime

A `Math::LinearSystem` is bound to the spaces that created it. A changed mesh,
space, block layout, or global size means a new system object. PETSc makes
this especially strict: matrices and vectors cannot be resized in place after
layout/assembly.

Warm starts are intentional. Reassembly zeroes the operator and right-hand
side as needed but preserves the solution vector where the backend contract
does so.

## Quadrature and exactness

Quadrature order is a numerical contract, not only a performance knob. A form
that claims exactness or a convergence rate must state the polynomial degree
being integrated and choose a rule that supports it. Curved mappings and
non-polynomial coefficients change that accounting.

When changing a fast path, compare against the generic quadrature path or a
manufactured solution. Identical-looking formulas can differ through mapping,
weights, or shape-function derivatives.

## Residual/tangent consistency

Every residual/tangent pair must satisfy a finite-difference consistency
check:

```text
J(x) w ~= (R(x + eps w) - R(x - eps w)) / (2 eps)
```

Run the check at P1 and at higher order if the implementation claims
order-genericity. A pair can be self-consistent and still physically wrong, so
also keep energy identities, manufactured solutions, or known-limit tests
when those are available.

## Constraints

Dirichlet conditions and identification conditions are structural constraints.
Assembly eliminates or expands them through the constraint map. Do not replace
them with penalty terms unless the mathematical model explicitly asks for a
penalty method.
