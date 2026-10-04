# P0 projection h-convergence

The real and complex scalar and two-component vector spaces are tested on all
seven `UniformGrid` cell geometries. For each value type, the discrete field
is obtained by solving the mass problem

$$
  (u_h,v_h)_\Omega=(f,v_h)_\Omega \qquad(v_h\in P_0).
$$

The exact field is affine and varies in every tested component. Its cellwise
means are the orthogonal projection onto the discontinuous constant space, so
the independently integrated L2 error decreases as $O(h)$. Three mesh levels
are checked. Complex tests use the modulus in the error norm and the library's
sesquilinear trial/test convention in the mass problem.

The same executable checks the global constant space P0g. Its exact
representing space contains every constant scalar or vector field, so the
assembled projection must reproduce four real/complex scalar/vector constants
to roundoff on each geometry. No h-rate is claimed for nonconstant P0g fields:
the global one-cell-independent space does not grow under mesh refinement.
