# Benchmark scope inventory


The middle column records available benchmark scopes, not exhaustive validated
or profiled coverage. Each scope has its own registrations and numerical
guards. The final column states the principal gaps; the presence of an
executable does not close them.

| Layer / owning module | Available measurement scopes | Principal remaining coverage |
| --- | --- | --- |
| Spatial and sparse algebra / Math | Scalar, vector and matrix contractions; mixed-scalar products; storage lifecycles; sparse products, dots and norms | Broader assembly-consumer attribution, mixed operators and distributed representations |
| Reference quadrature / QF | Construction, copying, traversal, repeated lookup and changes of rule identity/order | Process-first allocation, concurrent access and systematic order-dependent profiles |
| Reference elements / Variational, P1, H1 | Values and derivatives; scalar/vector/matrix basis production; prepared reference consumption; first binding and changing working sets | Complete controlled degree/order profiles, DOF functionals and component-expansion attribution |
| Transformations / Geometry | Pointwise and bulk Jacobians; construction, copy, move and rebinding; mapped coordinates, inverse and measure; mapped-quadrature lookup | Independent geometry-degree, curvature and embedding matrix; inverse-map stages |
| Mesh and space setup / Geometry, Variational, MPI | Mesh/connectivity construction, P1 setup and selected distributed construction | P0/P0g/H1 setup, restrictions, reconciliation and systematic degree scaling |
| Fields / Variational | Interpolation, mass projection, point evaluation, coefficient lookup, derivatives and reference contractions | Complete constant-space, matrix-range, high-degree and backend parity; norms and boundary restriction |
| Expressions and lifecycles / Variational | Named and generic forms, scalar operand order, basis scaling, Problem construction and assignment | Complete promotion/range matrix, cold expression construction and state-transition costs |
| Local integrators / Variational | Mass, weighted mass, directional derivative forms, gradient/Jacobian contractions and mixed Stokes bindings | Boundary/interface and nonlinear stage isolation; broader mixed degrees, embeddings and distributed consumers |
| Named mass kernels / Variational | Fresh-kernel and retained-reference preparation, coefficient replacement and component contractions | Component-expanded reference storage, other spaces and wider mapped/backend coverage |
| Local physical kernels / Solid | LinearElasticity first/repeated cell kernels, low-order and quadratic fields, affine and quadratic maps | Higher degrees, complex fields, coefficient variation and embedded formulations where applicable |
| Global assembly / Assembly | Complete operators/loads, sequential binding/kernel/scatter/completion stages, constrained scalar problems | Separate sparsity/setup costs; OpenMP merge, vector/mixed constraints and identification |
| Backend assembly / PETSc, MPI | Complete local/distributed physical, mass and constrained assembly consumers | Insertion/completion/communication attribution, ownership-sensitive costs and backend scaling |
| Linear solvers / Solver | Fixed-system CG and GMRES, diagonal setup/application, selected mixed direct solves | Factorization stages, broader preconditioners, distributed solves and setup/reuse separation |
| Complete workflows / owning modules | Interpolation/projection and assembly/solve consumer measurements | Matched setup/assembly/solve/postprocessing accounting across the applicable matrix |
| Supporting modules / Location, IO | Search and mesh-input benchmarks | Outside the assembly-first tranche; broader workflow attribution remains separate |
