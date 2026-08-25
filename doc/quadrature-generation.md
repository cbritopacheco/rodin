# Quadrature rule generation — state, gaps, and what was tried

Branch `feature/xg-wv`, 54 commits off `develop`. Everything below is generated
by code in this repository and checked against an independent moment oracle; no
coefficient was transcribed from a published table.

Two families are generated, and they are **not interchangeable**:

| family | construction | elements | file |
|---|---|---|---|
| Witherden–Vincent (WV) | fully symmetric, orbits of the element's symmetry group | all six | `src/Rodin/QF/WitherdenVincent.cpp` |
| Xiao–Gimbutas (XG) | node elimination, asymmetric | simplices only | `src/Rodin/QF/XiaoGimbutas.cpp` |

They publish **different counts for the same element and strength** — the XG
tetrahedron of strength 4 has 11 points and no symmetry, the WV one has 14 and
is fully symmetric. Comparing across families is meaningless and has already
caused one wasted day of work. Published counts are tabulated in
`tests/unit/Rodin/QF/PublishedCounts.h`.

---

## 1. What currently ships

### Witherden–Vincent (65 rules)

| element | strengths | counts |
|---|---|---|
| triangle | 1–20 | 1, 3, 6, 6, 7, 12, 15, 16, 19, 25, 28, 33, 37, 42, 49, 55, 60, 67, 73, 79 |
| quadrilateral | 1–15 | 1, 4, 4, 8, 8, 12, 12, 20, 20, 28, 28, 37, 37, 48, 48 |
| tetrahedron | 1–8 | 1, 4, 8, 14, 14, 24, 35, 46 |
| wedge (prism) | 1–7 | 1, 5, 8, 11, 16, 28, **39** |
| pyramid | 1–6 | 1, 5, 6, 10, 15, **23** |
| hexahedron | 1–9 | 1, 6, 6, 14, 14, 34, 34, 58, 58 |

All match the published count except: **pyramid strength 6 is 23 against a
published 24** (better — verified independently, see §6), and **wedge strength 7
is 39 against a published 35** (worse — the only shipped rule above published).

### Xiao–Gimbutas (29 rules)

| element | strengths | counts |
|---|---|---|
| triangle | 1–20 | 1, 3, **4**, 6, 7, **11**, **12**, 17, 19, **24**, **27**, **32**, **35**, **41**, **46**, **52**, **58**, **65**, **71**, 79 |
| tetrahedron | 1–9 | 1, 4, 6, 11, 15, 23, 31, 44, 57 |

Bold = **below** the published XG count: 13 of 20 triangle strengths. Strengths
1, 3, 5, 7, 9, 13, 15 sit exactly on the counting bound `⌈C(p+d,d)/(d+1)⌉`, so
those are **minimal** — no rule of that strength can have fewer points.

Two are above published: triangle strength 8 (17 vs 16), and tetrahedron
strength 5 (15 vs 14). The latter is **structural**: the 14-point rule attains
the counting bound exactly, which leaves the moment system square and its
solutions isolated, and node elimination works by moving along the null space of
an *under*determined system. Do not spend effort there without a different
method; the symmetric generator finds 14 directly, which is why WV tet 5 is 14.

---

## 2. Remaining gaps, in priority order

**The target is strength 20 on every element.** Against that, 66 of the 160
element-strengths are missing, not the nine rows below:

| element | have | missing |
|---|---|---|
| triangle (WV) | 20/20 | — |
| triangle (XG) | 20/20 | — |
| quadrilateral | 15/20 | 5 (16–20) |
| hexahedron | 9/20 | 11 (10–20) |
| tetrahedron (XG) | 9/20 | 11 (10–20) |
| tetrahedron (WV) | 8/20 | 12 (9–20) |
| wedge | 7/20 | 13 (8–20) |
| pyramid | 6/20 | 14 (7–20) |

The table below ranks only the strengths where a *published count* exists to
aim at, because that is where "did we meet it" has an answer. Beyond the
published range — strength 11+ in three dimensions, 16+ on the quadrilateral —
there is no target count, and a rule is judged by existing at all, by the
counting bound from below, and by the oracle. Do not read the absence of a
published number as the absence of a requirement; that framing is what made an
earlier version of this note stop ten strengths short on four elements.

Be aware of the scale involved. The counting bound at strength 20 in three
dimensions is **443 points**; the largest three-dimensional rule here is 58.
Nothing in this search has been exercised near that size, the per-solve cost at
that scale is untested, and the pyramid at strength 8 already enumerates 1.2M
decompositions. Strength 20 in three dimensions is not the same work with more
of it, and there is no evidence yet that this machinery reaches it.

Measured 2026-08-25. `cands@pub` is the number of decompositions at the
published point count; `≤6 orb` is how many of those use at most six orbits and
lie strictly inside the element.

| element | deg | published | now | bound | conditions | strata | cands@pub | ≤6 orb |
|---|---|---|---|---|---|---|---|---|
| wedge | 7 | 35 | 39 | 30 | 17 | 18 | 26 400 | 9 |
| tet | 9 | 59 | — | 55 | 18 | 11 | **708** | 2 |
| tet | 10 | 81 | — | 72 | 23 | 11 | **2 440** | 1 |
| hex | 10 | 90 | — | 72 | 16 | 14 | **1 136** | 18 |
| wedge | 8 | 46 | — | 42 | 24 | 18 | 149 232 | 4 |
| quad | 16 | 60 | — | 51 | 25 | 7 | 816 | 0 |
| quad | 17 | 60 | — | 57 | 25 | 7 | 816 | 0 |
| pyramid | 7 | 31 | — | 30 | 26 | 16 | 64 576 | 0 |
| pyramid | 8 | 47 | — | 42 | 35 | 16 | 1 234 752 | 0 |

**Start with wedge 7.** It is already proven reachable: the decomposition
`[2 3 6 6 6 12]` (orbit sizes; 17 unknowns against 17 conditions) solves at
residual **1.11e-15**. It has not been inlined only because the full search had
not yet reached it under the ordering in place at the time. The
fewest-orbits-first ordering committed in `a3f237174` should now find it; a run
of `RodinWitherdenVincent wedge 7 0 7` is the immediate next step.

**Then tet 9/10 and hex 10.** These have small candidate sets (708, 2 440,
1 136) — an exhaustive sweep is minutes, not hours. They were only ever
attempted under starved budgets or the round-robin bug (§5).

**Pyramid 7+ and quad 16+ are the hard ones.** Note the `≤6 orb` column is
**zero** for both: their orbits are small (pyramid 1/4/4/8, quad 1/4/4/8), so a
31- or 60-point rule needs 7+ orbits and the fewest-orbits heuristic gives no
purchase. Pyramid strength 8 has 1.2M decompositions. These likely need a
different idea, not more compute.

---

## 3. The stack

One generic pipeline; nothing is per-element except two small tables of data.

- **`SymmetryGroup.h`** — *derives* each element's affine symmetry group by
  trying maps that send vertices to vertices, then derives the orbit types as
  the fixed subspaces of subgroups (closed under intersection, deduplicated by
  group equivalence). Reproduces the published orbit lists for all six elements:
  triangle 1/3/6, quad 1/4/4/8, tet 1/4/6/12/24, prism 1/2/3/6/6/12, pyramid
  1/4/4/8, hex 1/6/8/12/24/24/48. Also derives facet inequalities, and the
  polytope each stratum's seed may occupy.
- **`CollapsedBasis.h`** — one orthogonal basis for *all seven* reference
  elements: the warped Jacobi tensor product (Karniadakis–Sherwin), driven by
  two small descriptor tables. Orthonormal to 1e-14 through degree 20 including
  the pyramid. Also supplies the exact positive collapsed Gauss–Jacobi rule used
  to measure norms, and a per-point tabulation for speed.
- **`SymmetricRuleGenerator.h`** — the WV search. Weights are eliminated (linear
  once points are fixed, solved by SVD each iteration); seeds are softmax
  combinations of the stratum polytope's vertices, hence interior for any
  parameter value and smooth everywhere; the Jacobian is analytic throughout,
  including through the elimination (Golub–Pereyra variable projection); modes
  that symmetry satisfies identically are pruned by probing.
- **`NodeElimination.h`** — the XG search. Unchanged this cycle; still uses
  `OrthogonalBasis.h` (simplex-only KD) and `SymmetricRuleSolver::MomentData`.

**Technical debt:** the old `SymmetricRuleSolver.h` still exists and is used by
`NodeElimination` for `MomentData`. Migrating `NodeElimination` to
`CollapsedBasis` would unify the stack and extend node elimination to
non-simplices, which may be the way into the pyramid/quad gaps.

---

## 4. Regenerating tables

```
RodinWitherdenVincent <tri|quad|tet|wedge|pyr|hex> [maxDegree] [seconds] [from] [jobs]
RodinXiaoGimbutas     <tri|tet>                    [maxDegree] [jobs]
```

Element is the only required argument; it names the array written. Defaults
reproduce the shipped tables: no deadline, whole machine, published range.

**`seconds` is the only argument that changes the answer.** It truncates the
search, so a strength that would have been solved comes back larger or not at
all. Zero (default) means no deadline. Use it only for unattended sweeps.

`jobs` changes how long a search takes, not what it returns — ordering and
seeding are fixed. It matters *through* the deadline, which is how several
strengths were wrongly recorded as unreachable.

Determinism is **per version of the search**: reordering candidates or changing
seed budgets can land on a different rule of the same size, equally exact. The
guarding test therefore asserts point counts, not coordinates.

---

## 5. What was tried — the log

### Worked

1. **Weight elimination** (WV's own method): weights are linear once points are
   fixed, so solve them by SVD each iteration instead of carrying them as search
   variables. ~⅓ fewer unknowns.
2. **Stick-breaking / softmax parameterisation** instead of clamping. A clamp
   has zero derivative once it binds, so any iterate starting outside the valid
   box had a collapsed Jacobian column and never moved again.
3. **Analytic Jacobian** (variable projection). A central difference caps the
   attainable residual near 1e-10 *however the step is tuned* — that cap was the
   accuracy of the rules, not the search.
4. **Orthogonal basis instead of Cholesky of the monomial Gram matrix.** The
   latter is Hilbert-like: marginal at degree 8, meaningless beyond.
5. **Positivity as part of the objective**, not a verdict afterwards. A
   configuration admits several roots and the bare moment residual routinely
   finds one with a negative weight.
6. **Pruning symmetry-trivial modes**, discovered by probing rather than derived
   per element.
7. **Molien's series** for the invariant dimension (Burnside on the action on
   polynomials, read off `1/det(I − tA)`). Exact and cheap.
8. **Ordering: fewest orbits before tightness of fit.** Published rules use a
   handful of large orbits; enumeration is dominated by many-small-orbit shapes
   that never close. Of 26 400 ways to make 35 points on a prism, 9 use ≤6
   orbits and one of those 9 is the rule.
9. **Caching basis norms** across decompositions, and **parallelising** over
   decompositions (3.4× wall-clock).
10. **Polish phase** after acceptance — a few extra Gauss–Newton steps take a
    rule from "clears the 1e-12 bar" to rounding.

### Failed, and why — do not repeat

1. **Treating the counting condition as necessary.** `unknowns ≥ invariant
   dimension` is **not** a necessary condition, and the rules that violate it are
   the economical ones. The cube's 6-point rule is one unknown short of two
   conditions; the pyramid's 23-point rule at strength 6 is *three* short of 20.
   Excluding such decompositions lost both. Cost: two regressions.
2. **Blaming the positivity floor.** Hypothesised that the 1e-4 relative weight
   floor was rejecting valid rules. Measured: every shipped rule sits at ≥3.5e-2
   of the mean weight, 350× above the floor. Hypothesis dead — do not revisit.
3. **Computing basis integrals with Grundmann–Möller.** Its weights alternate in
   sign; at these degrees cancellation put noise where exact zeros belong. The
   search then converged happily to residual 5e-14 **against the wrong
   equations**, producing rules wrong by 1e+02. Non-constant modes integrate to
   exactly zero by orthogonality — assert it, don't measure it.
4. **A `dJacobi` that could not express its own derivative.** The identity
   raises *both* Jacobi parameters, and the routine was fixed at β=0, so it
   silently returned the wrong polynomial. Gradients were wrong by O(1) until
   the general two-parameter recurrence went in.
5. **Round-robin that resampled instead of deepening.** Restarts are spread over
   candidates in rounds; each round was restarting the seed sequence from a
   different base, so four rounds of eight seeds never saw what one run of 32
   would have seen. The wedge's 28-point rule converges on the **21st** seed:
   8 restarts leave 3.7e-02, 32 reach 1.0e-15. This single bug suppressed
   *every* hard 3D strength and cost several hours; the strength was recorded as
   unreachable through an hour at full width, which is what a search that cannot
   deepen looks like from outside.
6. **Starving searches of cores while holding a deadline.** The hexahedron at
   strengths 6–9 was recorded as unreachable. An exhaustive sweep of every
   candidate there costs ~30 seconds; it had been given two threads and a
   15-minute deadline while three other sweeps competed. All four strengths were
   recovered in under 25 minutes once run properly. **Before concluding a
   strength is out of reach, give it the machine and no deadline.**
7. **Measuring the symmetric generator against the XG column.** Wasted
   considerable effort "fixing" a tetrahedron that was already at its own
   family's published count.
8. **Diagnostics that filtered out the answer.** Twice: a display filtered to
   `unknowns ≥ conditions` (hiding the underdetermined rules that work), and one
   limited to ≤5 orbits (hiding the 6-orbit rule that works). When a diagnostic
   says "nothing works here", check what it excluded.

---

## 6. How results are verified

Layered so a failure localises:

- **Rules as objects** — exactness against an independent monomial oracle,
  positive weights, interior points, invariance under the element's symmetry
  group, fuzz tests.
- **Counts** — every shipped rule asserted against its *own* family's published
  column, and against the counting bound from below (a count below the bound is
  not a discovery but a moment system that has stopped spanning). The three
  outstanding strengths are named in the tests with their reasons rather than
  silently excluded.
- **Machinery, independently of any rule** — basis orthonormality, spanning,
  gradient-vs-difference; group axioms, element preservation, orbit sizes vs
  published. Errors here do not surface as failures downstream, only as
  confident wrong answers.
- **Generator** — Molien dimension vs brute-force rank; defaults reproduce
  shipped counts.

**Better-than-published claims were verified outside the C++ entirely.** The
23-point pyramid rule was checked in Python against exact rational moments
derived by hand (`∫x^a y^b z^c = c!(a+b+2)!/((a+1)(b+1)(a+b+c+3)!)` on the
corner-apex reference pyramid): exact to 4.4e-15, least weight 3.4e-3, least
distance to boundary 1.5e-3. The XG triangle rules at strengths 3, 7, 13 were
checked in long double against exact rational moments: ~1e-15, strictly positive,
strictly interior, healthy margins.

Full rigour would need interval arithmetic; these are strong numerical
certificates, not proofs.

**Test count:** 140 quadrature tests, 288 geometry/variational/math, zero
failures.

---

## 7. Landmines

- The reference **pyramid is corner-apex**, not right — apex over vertex
  (0,0,0), not over the base centre. Its symmetries are therefore **shears, not
  orthogonal maps**, so any basis mapped through them must be re-orthonormalised
  before comparison. Missing this counted the two diagonal mirror planes as
  distinct orbit types.
- **Boundary orbits are required.** The only fully symmetric 6-point rule on a
  cube is its face centres — asking a seed `(t,½,½)` to integrate `x²` exactly
  forces `t` to 0 or 1. Excluding the boundary loses a rule that exists. They
  are admitted but tried after interior ones.
- `timeout` is **not available** on this macOS host; `gtimeout` or a wrapper is
  needed.
- Commits are **unsigned** — GPG cannot prompt in this shell. Run
  `git rebase --exec 'git commit --amend --no-edit -S' develop feature/xg-wv`.
- `clang-format` must be **18.1.8** (`/opt/local/libexec/llvm-18/bin`); v14
  silently ignores the repo config.
- Do not rebuild a target whose binary is currently running — the link fails.
