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
conjugates `rhs`. Every Rodin-defined `dot` member follows the same
convention, including the Frobenius contraction of dense matrices.
`Math::Vector`, `Math::FixedSizeVector`, `Math::Matrix` and
`Math::FixedSizeMatrix` own Eigen storage without adding data members.
Their inherited views, expressions and factory expressions retain Eigen
types and Eigen member-dot semantics. Use `Math::dot` for those operands.
The matrix classes derive from `Eigen::PlainObjectBase` with the Rodin
owner as the CRTP argument. Deriving from a concrete `Eigen::Matrix` would
retain Eigen's native derived type and can invalidate references inside
explicit Eigen decomposition templates. Use `Eigen::MatrixBase` interfaces
or `Eigen::Ref` to borrow a Rodin matrix in backend code. Constructing a native
Eigen matrix copies its value; moving matching native dynamic storage into a
Rodin matrix transfers the allocation.
A complex-valued form has the discrete actions

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
`conj(c) * conj(u) * v`. The local backend evaluates the action as
`Math::dot(A * u, v)` and PETSc as `VecDot(A * u, v)`. Both conjugate the
second operand. Raw Eigen member dot conjugates its first operand and would
instead require `v.dot(A * u)`. The corresponding linear-form actions are
`Math::dot(b, v)` and `VecDot(b, v)`. Do not infer a backend's convention
from a member method's spelling alone.

Sesquilinearity is not Hermitian symmetry. Whether `A == A*` is a separate
property of the particular weak form and coefficients.

`LinearElasticityIntegral` admits Lamé parameter ranges whose multiplication
preserves the field scalar type. Real parameters may therefore act on complex
fields, but complex values must not be narrowed into real storage. Geometry
inverses remain real; complex basis Jacobians use the existing promoted matrix
product, while the established real in-place arithmetic is retained.

For coincident trial/test spaces, triangular reconstruction is valid when both
Lamé parameter ranges are real. Complex-valued material ranges do not establish
$A=A^*$, even when the same space is used twice, so both matrix triangles are
evaluated. For the currently real-valued basis functions, complex material
parameters may yield $A^T=A$ without Hermitian symmetry. P1 and generic
integrators follow the same rule; it is independent of the assembly backend.

## LinearSystem lifetime

A `Math::LinearSystem` is bound to the spaces that created it. A changed mesh,
space, block layout, or global size means a new system object. PETSc makes
this especially strict: matrices and vectors cannot be resized in place after
layout/assembly.

Warm starts are intentional. Reassembly zeroes the operator and right-hand
side as needed but preserves the solution vector where the backend contract
does so.

## Physical differential dimensions

Let a reference cell have dimension $d$, its physical embedding have dimension
$s$, and a vector field have $m$ components. The geometric inverse has shape
$d\times s$. Physical scalar gradients and vector Jacobians therefore have
shapes $s$ and $m\times s$, respectively:

$$
\nabla_K u=(J_K^{\dagger})^T\widehat\nabla u,
\qquad J_K(u)=\widehat J(u)J_K^{\dagger}.
$$

For an embedded cell these are tangential derivatives expressed in ambient
coordinates. Reference-coordinate contractions still have $d$ terms; derivative
outputs retain all $s$ ambient entries. Neither the field component count nor
the integration entity's topological dimension determines the output width.
Constant-space derivatives obey the same shape contract even though all entries
are zero. On a point embedded in $\mathbb{R}^s$, the empty reference contraction
produces a correctly shaped ambient zero, not an empty physical differential.
Direct point evaluation, indexed integration points, shape-function evaluation
and specialized local integrators must agree on these dimensions.

Global-constant derivative nodes borrow immutable operands and update only their
own cursor and zero storage. Their derivative does not require mutating the
operand's evaluation cursor.

## Quadrature and exactness

Quadrature order is a numerical contract, not only a performance knob. A form
that claims exactness or a convergence rate must state the polynomial degree
being integrated and choose a rule that supports it. Curved mappings and
non-polynomial coefficients change that accounting.

When changing a fast path, compare against the generic quadrature path or a
manufactured solution. Identical-looking formulas can differ through mapping,
weights, or shape-function derivatives.

### Family order, degree and geometry

Concrete `P0`, `P0g`, `P1` and `H1` elements report their family index through
`getOrder()`. `FiniteElementBase::getDegree()` instead reports an optional
upper bound on total polynomial degree of reference basis components.
Tensor products therefore have different orders and degrees. The overload
`getDegree({directions...})` bounds the joint degree in selected reference
coordinates, treating unselected coordinates as coefficients. On a wedge,
`{0, 1}` preserves the triangular total-degree bound and `{2}` the segment
bound. Invalid directions and unknown bounds return an empty optional;
repeated directions are ignored. Rational pyramid bases have a family order
but no established total polynomial degree bound; their horizontal coordinates
have known bounds when the vertical coordinate is treated as a coefficient.
The shared directional Lagrange-space calculation belongs to
`LagrangeFiniteElementBase`, which does not impose nodal data or a particular
basis representation. P0 and P0g validate directions and report their constant
degree directly; unrelated finite element families do not inherit this policy.

Expressions retain `getOrder(polytope)` as an integration recommendation,
including documented nonpolynomial estimates. Zero additionally certifies an
elementwise constant expression; constant-value caches rely on this contract.
Do not substitute a family index for an existing total-degree estimate without
examining the operator, geometry and selected rule. Generic rules can consume
`getDegree()` without requiring a family index. Quadrature selection belongs to
expressions and the applicable rules; native P1/H1 rules retain their documented
rational-pyramid estimates. If neither the expression nor the rule establishes
an order, require an explicit integration override. Do not replace an unknown
estimate with a universal numerical fallback.

Transformations expose total and directional `getDegree()`, `isFirstOrder()`
and `getBounds()`.
Native first-order product and rational pyramid maps need not be affine.
`Bounds` owns a conservative convex enclosure with outward support bounds;
its pyramid construction uses certified collapsed-reference-space containment,
never a polynomial inference from samples. Unsupported maps remain unpruned.

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

### Targeted assembly with prescribed values

For a square system with matching trial/test constraint indices, let
\(P\) project onto unconstrained DOFs and let \(g\) contain the prescribed
values, with zero entries on unconstrained DOFs. Full assembly produces

$$
A_D = PAP + (I-P), \qquad b_D = P(b-Ag)+g.
$$

An LHS-only update produces the same row-and-column-eliminated operator
\(A_D\), leaving the stored right-hand side untouched. An RHS-only update
produces \(b_D\), leaving the stored operator untouched. Its lifting uses
the current unconstrained bilinear form, including local, nonlocal and
preassembled contributions, rather than the entries of a previously
eliminated or externally scaled stored matrix. Successive target updates
therefore reconstruct the same constrained system as full assembly.


### Affine identification over real and complex fields

For an affine expansion \(x=T\hat{x}+d\), assembly distributes test rows
through the conjugated expansion and trial columns through the original
expansion. The constrained form and load are \(T^{*}AT\) and
\(T^{*}(b-Ad)\), followed by explicit reconstruction rows for retained slave
unknowns. Over real fields the adjoint reduces to the transpose. Reconstruction
coefficients retain their original phase.

PETSc targeted assembly rejects identification constraints before resetting
any system output. For MPI assembly this rejection is collective: every rank
rejects when any rank has identification rows, including ranks with empty
shards. The collective belongs to the explicit partial-assembly operation;
full assembly and local accessors do not gain this validation collective.
