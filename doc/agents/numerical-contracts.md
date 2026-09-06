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

