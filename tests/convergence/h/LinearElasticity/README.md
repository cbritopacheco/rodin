# LinearElasticity h-convergence

This suite validates the conforming displacement formulation of isotropic
linear elasticity,

$$
  -\nabla\cdot\sigma(u)=f, \qquad
  \sigma(u)=\lambda(\nabla\cdot u)I+2\mu\varepsilon(u), \qquad
  \varepsilon(u)=\tfrac12(\nabla u+\nabla u^T).
$$

For a smooth displacement represented by degree-`K` conforming elements,
Cea's lemma, Korn's inequality, interpolation estimates, and the usual duality
argument predict

$$
  \lVert u-u_h\rVert_{H^1}=O(h^K), \qquad
  \lVert u-u_h\rVert_{L^2}=O(h^{K+1}).
$$

The error oracle integrates the Euclidean displacement error and the Frobenius
norm of the displacement-Jacobian error independently of the assembled form.
The tests cover:

1. exact P1 reproduction of an affine displacement and constant strain;
2. optimal P1 and P2 rates for a smooth exponential displacement with full
   Dirichlet data;
3. the same rates with displacement imposed on `x_0=0` and manufactured
   traction `sigma(u)n` on the complementary boundary;
4. P2 convergence of the transverse shear wave
   $u=(\sin(\pi x_1),0,\ldots,0)$ with `lambda=10^4` and `mu=1`. Both the
   exact field and its element-wise nodal interpolants have zero divergence,
   so this checks the nearly incompressible parameter regime without allowing
   volumetric locking to obscure the approximation rate.

Segment, triangle, quadrilateral, tetrahedron, pyramid, hexahedron, and wedge
are covered wherever the mathematical case applies. The nearly incompressible
divergence-free case excludes Segment because a nonzero divergence-free
displacement does not exist in one dimension. Boundary normals on the segment
use the exact endpoint orientation; higher-dimensional tractions use Rodin's
`BoundaryNormal`.

The mixed P1 tetrahedron study starts one refinement level later than the
other geometries. Its coarsest unit-box grid is visibly pre-asymptotic in the
L2 dual estimate at the displacement/traction interface; the documented
levels ensure that both measured intervals test the asymptotic regime.
