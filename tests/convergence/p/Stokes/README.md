# Stokes degree-refinement verification

For $\Omega=(0,1)^d$, $d\in\{2,3\}$, viscosity $\nu=1$, and coordinates
$x_0,\ldots,x_{d-1}$, the manufactured Stokes problem is

$$
-\nu\Delta u+\nabla p=f,\qquad \nabla\cdot u=0,\qquad
u|_{\partial\Omega}=g,\qquad \int_\Omega p\,dx=0.
$$

The analytic fields and source, supplied by `StokesData`, are

$$
u=\sin(\pi x_1)e_0,\qquad p=\cos(\pi x_0),\qquad
f=\bigl(\pi^2\sin(\pi x_1)-\pi\sin(\pi x_0)\bigr)e_0,
$$

where $e_0$ is the first coordinate unit vector and $g$ is the full exact
velocity trace. The only nonzero Jacobian entry is
$(\nabla u)_{01}=\pi\cos(\pi x_1)$; the pressure gradient is
$\nabla p=-\pi\sin(\pi x_0)e_0$. These references are evaluated in physical
coordinates independently of the discrete solutions. Nonpolynomial data
prevents finite-degree exact reproduction from replacing the refinement study.

## Mixed spaces and mathematical scope

The velocity degree $k$ and pressure degree $k-1$ increase together through
$(k,k-1)=(2,1),(3,2),(4,3)$ on a fixed mesh with `n=3` points per axis,
$h=1/2$. The local scalar/vector H1 families supply the cell-appropriate
spaces; a real $P_0^g$ multiplier fixes the pressure mean. The degree offset
is the [Taylor–Hood construction](https://defelement.org/elements/taylor-hood.html)
on simplices, with tensor-product analogues on quadrilateral/hexahedral cells.
The same degree offset is tested on Rodin's pyramid and wedge spaces; it
is not, by itself, a stability theorem for those geometries.

For homogeneous velocity variations $V_{h,k}^0$ and zero-mean pressures
$Q_{h,k-1}^0$, mixed well-posedness requires a positive discrete inf-sup
constant

$$
\beta_{h,k}=\inf_{0\ne q_h\in Q_{h,k-1}^0}
\sup_{0\ne v_h\in V_{h,k}^0}
\frac{|(q_h,\nabla\cdot v_h)_{L^2}|}
{\lVert q_h\rVert_{L^2}\lVert\nabla v_h\rVert_{L^2}}>0.
$$

Under the relevant stability and regularity hypotheses, approximation
improvement is expected as both degrees grow. The tests certify specified
finite workloads and field errors; they do not establish a uniform lower
bound on $\beta_{h,k}$ or an asymptotic exponential theorem as $k\to\infty$.
In particular, simplex results are not silently transferred to pyramids.

The fixed mesh starts at `n=3`, not `n=2`. An exact DOF-count regression
assembles homogeneous velocity constraints on a single quadrilateral or
hexahedron. For the degree-two/degree-one tensor-product pair,

$$
\dim V_{h,2}^0=d<2^d-1=\dim Q_{h,1}^0,\qquad d\in\{2,3\}.
$$

The discrete divergence cannot have full pressure rank, so a nonconstant
pressure null mode remains even after fixing the mean. This is a logical
rank obstruction, checked through actual constrained indices without a
floating-point rank tolerance. Avoiding that coarse mesh does not, by
itself, prove a uniform inf-sup bound on the larger meshes.

`StokesProblem` is shared by native h, p, and hp studies. Each solve creates
fresh spaces and a fresh saddle-point system. In the form language it states

$$
\begin{aligned}
\nu(\nabla u_h,\nabla v_h)-(p_h,\nabla\cdot v_h)&=(f,v_h),\\
(\nabla\cdot u_h,q_h)+\lambda_h(1,q_h)&=0,\\
(p_h,1)&=0.
\end{aligned}
$$

The multiplier does not make the velocity pointwise divergence-free. Strong
L2 divergence is checked for exactly representable patches, not required to
decrease monotonically in the analytic study.

## Quantities and acceptance

Independent physical-cell quadrature measures four errors:

$$
E_{u,0}=\lVert u-u_h\rVert_{L^2},\quad
E_{u,1}=\lVert\nabla u-\nabla u_h\rVert_{L^2},\quad
E_{p,0}=\lVert p-p_h\rVert_{L^2},\quad
E_{p,1}=\lVert\nabla p-\nabla p_h\rVert_{L^2}.
$$

Every error must be finite and positive, strictly decrease on both degree
intervals, and satisfy $\log(E_{k-1}/E_k)>0.1$. Since the velocity-degree
increment is one, this is the reported finite-degree logarithmic decay.
All four quantities are asserted separately; accurate velocity cannot conceal
inaccurate pressure. Three pairs provide two measured intervals.

SparseLU is used for the saddle-point system; no CG/SPD assumption is made.
Factorization/solve success and an independently recomputed coefficient
residual $\lVert Ax-b\rVert_2/\max(1,\lVert b\rVert_2)<10^{-11}$ are
required. Assembly and error quadrature use order 16. Pressure mean is
integrated at order 18 and must have absolute magnitude below $10^{-10}$.
The P4/P3 sensitivity case on the same `n=3` mesh raises assembly/error
quadrature to 18; each error must change by less than $10^{-6}$ relative to
baseline. SparseLU has no iterative residual-tolerance parameter to tighten.

## Polynomial and physical negative controls

For each pair, a representable patch excites its highest field degrees:

$$
u=x_1^m e_0,\qquad p=x_0^{m-1}-\frac1m,\qquad
f=\bigl(-m(m-1)x_1^{m-2}+(m-1)x_0^{m-2}\bigr)e_0,
\qquad m\in\{2,3,4\}.
$$

The respective pairs use $m=k$. All four field errors and strong L2
divergence must be below $10^{-9}$.

A separate control retains the quadratic reference $u=x_1^2e_0$,
$p=x_0-1/2$, source $f=-e_0$, and exact trace, but changes viscosity to
$\nu=2$. The wrong problem has the same velocity and pressure
$p_{\rm wrong}=3x_0-3/2$. Consequently,

$$
\lVert p_{\rm wrong}-p\rVert_{L^2}=\frac1{\sqrt3},\qquad
\lVert\nabla(p_{\rm wrong}-p)\rVert_{L^2}=2.
$$

At every tested pair, velocity errors must remain below $10^{-9}$ while
pressure errors must exceed $0.1$ and $1$, respectively. The pressure oracle
therefore demonstrably rejects the incorrect physics rather than relying
only on algebraic solver convergence.

## Execution scope

Triangle, quadrilateral, tetrahedron, pyramid, hexahedron, and wedge meshes
are tested with native Eigen sequential/OpenMP assembly. Point and segment
are excluded from this non-degenerate incompressible Stokes workload.
PETSc/MPI p studies, curved maps, and arbitrary higher degree are separate
gaps, not implied by native success. The existing PETSc h suite has its own
backend/rank matrix.

Each geometry entry has labels `convergence;slow` and a 600-second limit.
Run `ctest --test-dir build/tests -R RodinConvergencePStokes --output-on-failure`.
