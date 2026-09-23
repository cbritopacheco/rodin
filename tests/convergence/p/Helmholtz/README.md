# Complex Helmholtz p-convergence

The complex plane wave $u(x)=\exp(i\sum_j x_j)$ solves
$-\Delta u-\tfrac14u=(d-\tfrac14)u$ on the unit box, with its exact trace
prescribed on the boundary. The bilinear form contains both complex stiffness
and a negative mass term; the wavenumber is low enough that the Dirichlet
problem remains coercive on this domain.

On a fixed two-points-per-axis `UniformGrid`, this suite measures independent
L2 and H1-seminorm errors at polynomial degrees 1 through 4. Analytic
solutions admit exponential best-approximation decay in degree, subject to
the finite degree range and quadrature/solver errors. The test checks every
successive reduction and its observed exponential rate, on all seven cell
geometries. It complements the same equation's P1/P2 h-convergence suite.

As a negative control, removing the Helmholtz mass term while keeping the
manufactured forcing causes the degree-2-to-3 and degree-3-to-4 rate checks
to fail: the computed solution approaches a different PDE solution rather
than the prescribed plane wave. The term was restored after this check.
