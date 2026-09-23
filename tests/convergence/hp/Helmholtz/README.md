# Combined hp refinement for complex Helmholtz

The analytic plane wave $u(x)=\exp(i\sum_j x_j)$ satisfies
$-\Delta u-\tfrac14u=(d-\tfrac14)u$ with its exact complex trace imposed
on the boundary. The test uses the same manufactured equation as the h and p
Helmholtz suites, allowing the three refinement mechanisms to be compared.

The mesh widths halve while the degree increases: three all-geometry runs use
$(h,p)=(1,1),(1/2,2),(1/4,3)$ on `UniformGrid` meshes. Each run solves the
complex weak form and measures L2 and H1-seminorm errors by independent
quadrature. Both errors must strictly decrease at every step; their
log-ratios with respect to $h$ must exceed the P1 h-only benchmark powers of
two and one, respectively. These ratios measure the selected combined path,
not a universal fixed-order h exponent or a pure p rate.
