# PETSc variable-conductivity degree refinement

On $\Omega=(0,1)^d$, let $s(x)=\sum_jx_j$ and
$\gamma(x)=1+s(x)>0$. The real field satisfies

$$
-\nabla\cdot(\gamma\nabla u)=f,\qquad
u|_{\partial\Omega}=u_*,\qquad u_*=e^s,\qquad
f=-d(1+\gamma)e^s.
$$

The coefficient-gradient contribution is retained in the source. The
fixed $n=2$ grid, degrees $1\to2\to3\to4$, independent field norms,
every-interval decay floor $\log(E_{p-1}/E_p)>0.25$, and backend/rank
matrix match the native conductivity degree study. Analytic data and
positive conductivity support degree decay under the usual approximation
and integration hypotheses; the finite observed floor is not a theorem
with a universal exponential constant.

The shared problem/study classes, quadratic degree-threshold patch,
highest-degree reproduction, separate assembly/norm/solver sensitivity
controls and independent residual budget are specified in the
[Poisson counterpart](../PETScPoisson/README.md). They are reused, not
independently implemented for each coefficient or refinement axis.

The resolved wrong-coefficient control instead uses $u_*=1+s$ on $n=3$,
for which $f=-d$. The correct variable-coefficient operator first
reproduces this affine field within $10^{-9}$ in both norms. Replacing
$\gamma$ by $1$, while retaining $f$ and the trace, must produce
$E_{L^2}>10^{-3}$ and $E_{H^1}>10^{-2}$. This checks the operator's
variable coefficient, rather than accepting two consistently altered
manufactured problems. All geometry/rank registrations share the same
slow labels, safety timeouts and pyramid scheduling lock as Poisson.
