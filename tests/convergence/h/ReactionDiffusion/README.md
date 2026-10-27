# Coupled reaction–diffusion h-convergence

Two fields satisfy $-\Delta u+u+\alpha w=f_1$ and
$-\Delta w+w+\alpha u=f_2$, with $\alpha=0.2$. The exact fields are
$u=\prod_i\sin(\pi x_i)$ and $w=\prod_i\cos(\pi x_i)$; their sources
are $f_1=(d\pi^2+1)u+\alpha w$ and
$f_2=(d\pi^2+1)w+\alpha u$. Exact Dirichlet data are imposed on every
boundary, including the nonzero cosine field. The coupling matrix is
positive definite because $|\alpha|\lt 1$.

Each field is checked independently in L2 and H1 seminorm on three grids.
Conforming P1 and P2 should give rates $(2,1)$ and $(3,2)$, respectively.
All seven UniformGrid cell geometries are covered. This detects coupling,
field ordering, and nonhomogeneous boundary regressions that scalar Poisson
cannot expose.
