# Variational formulation — theory ↔ code

How the mathematics of weak forms becomes Rodin source. Read with
[../variational.md](../variational.md) (the inventory) at hand.

## From strong to weak form

A boundary-value problem

$$-\,\mathrm{div}(A\nabla u) = f \ \text{in}\ \Omega, \qquad u = g \ \text{on}\ \Gamma_D, \qquad A\nabla u\cdot n = h \ \text{on}\ \Gamma_N$$

is multiplied by a test function $v$ vanishing on $\Gamma_D$ and integrated
by parts:

$$\text{find } u \in V_g:\quad a(u,v) := \int_\Omega A\nabla u\cdot\nabla v
\;=\; \int_\Omega f v + \int_{\Gamma_N} h v =: \ell(v) \quad \forall v\in V_0.$$

In Rodin this *is* the program text:

```cpp
Problem poisson(u, v);
poisson = Integral(A * Grad(u), Grad(v))
        - Integral(f, v)
        - BoundaryIntegral(h, v).over(GammaN)
        + DirichletBC(u, g).on(GammaD);
```

Correspondences to keep straight:

- The statement is $a(u,v) - \ell(v) = 0$, so **load terms carry a minus
  sign** — the `-` in front of `Integral(f, v)` is mathematics, not style.
- **Essential (Dirichlet) conditions** restrict the space ($u \in V_g$) and
  are imposed on DOFs by elimination (`DirichletBC` →
  `Assembly::ConstraintMap` → `Math::LinearSystem` elimination). They never
  appear as integrals.
- **Natural (Neumann/Robin) conditions** arise from integration by parts
  and *are* integrals (`BoundaryIntegral`). Getting this backwards
  (penalizing Dirichlet data, or "setting" Neumann DOFs) is the classic
  formulation error.

## Function spaces and conformity

The continuous spaces are Sobolev spaces: $H^1(\Omega)$ (square-integrable
value and gradient), $L^2(\Omega)$, and — planned in the design though not
yet realized as spaces — $H(\mathrm{div})$ and $H(\mathrm{curl})$. A
discrete space is *conforming* when $V_h \subset V$:

- `H1<K>` / `P1` are $H^1$-conforming: continuous across facets, so
  `Grad` is square-integrable and `Integral(Grad(u), Grad(v))` is
  well-defined.
- `P0` (and DG use of facet integrals) live in $L^2$: no continuity, so
  derivatives exist only element-wise, and coupling between elements must
  be written explicitly through facet terms (below).

The trace theorem ($H^1$ functions have boundary values in
$H^{1/2}(\partial\Omega)$) is what makes `traceOf`, `TraceOperator` and
boundary integrals meaningful; `UndeterminedTraceDomainException` is the
code's complaint when a function is evaluated on an interior facet without
saying *which side's* trace is meant.

## Well-posedness (what you're implicitly assuming)

Lax–Milgram: if $a$ is bounded and coercive on $V_0$, the weak problem has
a unique solution. Practical consequences in code:

- Coercivity failing ⇒ the assembled matrix is singular/indefinite. Usual
  causes: no Dirichlet condition anywhere (pure Neumann ⇒ constant
  nullspace — fix the constant, e.g. with a `P0g` multiplier), or a sign
  error in a term.
- Symmetric coercive $a$ ⇒ SPD matrix ⇒ `CG`/`SimplicialLDLT` are valid.
  Nonsymmetric terms (advection, some couplings) ⇒ `BiCGSTAB`/`GMRES`/LU.
  Solver choice is a *theorem*, not a preference.

## The Galerkin method is a projection

Replacing $V$ by $V_h = \mathrm{span}\{\phi_i\}$ turns the weak problem
into $\sum_j a(\phi_j,\phi_i)\, U_j = \ell(\phi_i)$ — i.e. the stiffness
matrix $A_{ij} = a(\phi_j, \phi_i)$ and load $b_i = \ell(\phi_i)$. Céa's
lemma says the Galerkin solution is (up to a constant) the best
approximation in $V_h$; convergence rates then come from interpolation
theory: for degree-$K$ elements and a smooth solution,
$\|u-u_h\|_{H^1} = O(h^K)$ and $\|u-u_h\|_{L^2} = O(h^{K+1})$.
**tests/manufactured asserts exactly these orders** — a failing
manufactured test usually means a quadrature order, mapping, or conformity
bug, not a "tolerance issue".

In Rodin's type system the grading of the algebra mirrors the tensor
structure of forms (see philosophy.md): a grade-1 integrand is an element
of $V_h'$ (assembles to a vector over test DOFs), a grade-2 integrand an
element of $(U_h \otimes V_h)'$ (assembles to a matrix). `Integrator::Type
{Linear, Bilinear}` is this distinction at runtime; the
Trial/Test `ShapeFunctionSpaceType` tags are it at compile time.

## Facet terms: jumps, averages, DG

For discontinuous fields, element-boundary information enters through the
facet trace operators
$[\![w]\!] = w^+ - w^-$ (`Jump`) and $\{w\} = \tfrac12(w^+ + w^-)$
(`Average`), with `FaceNormal` fixing the orientation convention. Interior
penalty DG for $-\Delta u = f$ reads schematically

$$\sum_K \int_K \nabla u\cdot\nabla v
- \int_{\mathcal F_i} \{\nabla u\}\cdot n\, [\![v]\!]
- \int_{\mathcal F_i} \{\nabla v\}\cdot n\, [\![u]\!]
+ \frac{\eta}{h}\int_{\mathcal F_i} [\![u]\!]\,[\![v]\!],$$

written with `Integral` over cells plus `InterfaceIntegral` (interior
facets) / `BoundaryIntegral` (exterior facets); `FaceIntegral` ranges over
both. The penalty scaling $\eta/h$ is what makes the form coercive —
dropping it is the standard DG instability. See examples/DG.

## Constraints: elimination, identification, multipliers

Three mathematically distinct ways to impose $Bu = g$, all present:

1. **Elimination** (value `DirichletBC`): constrained DOFs are removed from
   the system; exact, no parameters. Preferred whenever the constraint is
   nodal.
2. **Identification** (`DirichletBC(u, A(v))`, `PeriodicBC`): boundary DOFs
   of one field are declared *equal to a linear expression* in another
   still-unknown field; the `ConstraintMap` performs the substitution at
   assembly. This is exact algebraic gluing (mortar-like, periodic,
   rigid-coupling conditions).
3. **Lagrange multipliers / global constraints**: an extra unknown enforces
   the constraint weakly — the natural home of the `P0g` global-constant
   space (e.g. zero-mean pressure, prescribed total volume). Produces
   saddle-point systems: SPD solvers no longer apply (use MINRES/LU), and
   inf-sup (LBB) compatibility between the spaces becomes the relevant
   well-posedness condition.

Penalty formulations (adding $\frac{1}{\epsilon}\int (u-g)v$) trade
exactness for conditioning and are *not* the house style for essential
constraints — though smooth penalties are exactly right for *soft*
constraints like surface fitting (conventions.md).

## Nonlinear problems

A nonlinear weak form $R(u; v) = 0$ is solved by Newton: at iterate $u^k$,
solve the linearized problem $R'(u^k)(\delta u, v) = -R(u^k; v)$ — a
bilinear form in $(\delta u, v)$ called the (consistent) tangent. In code:
a residual `LinearFormIntegrator` in $v$ paired with a tangent
`BilinearFormIntegrator` in $(\delta u, v)$, driven by
`Solver::NewtonSolver` (see [../solvers-assembly.md](../solvers-assembly.md)
for its damping/line-search limitations). The FD-consistency test pattern
(testing.md) is precisely the check that the tangent is the derivative of
the residual — Newton's quadratic convergence is the *symptom*; FD
consistency is the *property*.

## Integral equations

`Potential` applies a kernel operator
$(\mathcal K u)(x) = \int_\Gamma k(x,y)\, u(y)\, dy$ — the single-layer /
Newtonian potential machinery of boundary integral equations
(examples/IntegralEquations). Unlike differential terms these couple all
DOFs (dense blocks); they exist in the language as first-class nodes but
do not enjoy the sparsity assumptions of the standard assembly path.
