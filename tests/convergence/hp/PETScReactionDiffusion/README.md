# PETSc coupled reaction–diffusion combined refinement

The equation, two-field manufactured data, shared workload, residual
checks and componentwise controls are those of the
[degree-refinement suite](../../p/PETScReactionDiffusion/README.md).

The three-level path is

$$
(n,p)=(2,1)\to(3,2)\to(5,3),\qquad h=1/(n-1).
$$

For each field $j\in\lbrace 0,1\rbrace $ and each interval, the finite, positive
$L^2$ error $E_{i,j}^{(0)}$ and $H^1$-seminorm error $E_{i,j}^{(1)}$
must strictly decrease and satisfy

$$
r_{i,j}^{(s)}=\frac{\log(E_{i-1,j}^{(s)}/E_{i,j}^{(s)})}
 {\log(h_{i-1}/h_i)},\qquad
r_{i,j}^{(0)}>1.9,\quad r_{i,j}^{(1)}>0.9.
$$

These are effective slopes along this combined mesh/degree path, not
fixed-degree h orders or a universal hp theorem. No field or interval
can be hidden by an aggregate norm or endpoint-only comparison.

Degree three on $n=3$ supplies the quadratic patch, resolved affine
wrong-coupling control, and three separate quadrature/solver checks.
Their absolute and relative budgets are unchanged from the p suite.
All seven positive-dimensional cell families have local and MPI
rank-one-through-four registrations with sequential and OpenMP assembly.
Shared slow labels, timeouts and resource locks affect scheduling only,
not mathematical acceptance.
