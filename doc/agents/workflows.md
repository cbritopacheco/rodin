# Extension workflows

Use this page when the task is "add a new X" or "extend X to a new case".
The philosophy page explains the shape of the library; this page is the
operational checklist.

## New form-language node

1. Add one header for the algebraic concept, named for the mathematical
   operation (`Trace`, `Jump`, `Potential`, etc.).
2. Add `FormLanguage::Traits` specializations in the same header, above the
   class definitions.
3. Specialize the existing concept name for each grading/range case instead
   of creating sibling class names.
4. Deep-copy operands on construction and keep evaluation on the CRTP/static
   path. Do not add virtual calls to quadrature-point evaluation.
5. Add CTAD/deduction helpers so user code reads like the formula.
6. Add Doxygen math for the formula and a specialization table when the node
   has multiple supported specializations.
7. Test the expression in the smallest form that exercises its grading:
   function, linear-form integrand, bilinear-form integrand, or DG/facet term.

## New finite element space

1. Define the element first: basis, nodes/functionals, order, and range.
2. Define how DOFs distribute over the incidence complex by dimension. Never
   assume a `node * vdim + component` layout in users of the space.
3. Implement pullback/pushforward and point evaluation through
   `Geometry::Point` and the element basis.
4. Add `FormLanguage::Traits` for the space and its element.
5. Add FES-specific operator specializations (`Grad`, `Div`, `Jacobian`,
   `QuadratureRule`) only where the space can support them mathematically.
6. Add `GridFunction` coverage: assignment/interpolation, evaluation at a
   point, vector-valued behavior if applicable, and I/O if the space is meant
   to be saved.
7. Verify with manufactured convergence when the space claims approximation
   order, not just with compile tests.

## New solver

1. Start from `Solver::LinearSolverBase<LinearSystem>` unless the solver is
   truly nonlinear.
2. Add `FormLanguage::Traits<SolverType>` exposing `LinearSystemType`.
3. Keep the user surface `Solver(problem).setX(...).solve()` and preserve
   chainable setters.
4. Specialize on the linear-system type: sparse, dense, or PETSc-backed.
5. Document the matrix assumptions explicitly: SPD, symmetric indefinite,
   nonsymmetric, rectangular, direct factorization, or Krylov iteration.
6. Add a specialization table on the primary class page.
7. Test both the raw linear-system path and the `Problem`/CTAD path when both
   are supported.

## New assembly behavior

1. Identify the operand: bilinear form, linear form, Dirichlet value,
   identification constraint, or complete problem.
2. Add or extend the typed input object, preserving both local integrator
   lists and global integrator lists.
3. Implement `AssemblyBase` first, then concrete iteration backends:
   `Sequential`, `OpenMP`, MPI, and PETSc mirrors as applicable.
4. Preserve constraint semantics: Dirichlet and identification constraints
   are structural eliminations, not penalties.
5. Preserve `Math::LinearSystem` lifetime rules: a changed mesh or space means
   a new system object.
6. Add tests that cover the backend matrix the feature claims to support.

## New IO support

1. Extend the existing loader/printer family:
   `MeshLoader`, `MeshPrinter`, `GridFunctionLoader`, or
   `GridFunctionPrinter`.
2. Specialize by `IO::FileFormat`, context, finite element space, and data
   backend. Do not add a parallel naming hierarchy.
3. Decide whether the format is mesh-only, field-only, local, MPI, PETSc, or
   time-series output.
4. Round-trip the smallest representative mesh/field and verify attributes,
   geometry, connectivity requirements, and coefficient data.
5. Add the specialization table rows to the family page.

## New Solid term

1. Write the residual and the consistent tangent together.
2. State the sign convention in the docs. `Problem` moves linear forms to the
   right-hand side, so Newton systems should assemble `K du = -R`.
3. Keep internal variables as finite element unknowns when they participate
   in the coupled solve.
4. Add finite-difference consistency tests for the residual/tangent pair.
5. Run the test at P1 and at higher order if the term claims order-genericity.
6. For hyperelasticity, state the kinematics, stress measure, energy, and
   tangent being implemented.

