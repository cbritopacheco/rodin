/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @example Stokes problem with PETSc (Taylor–Hood P2–P1 + mean-pressure constraint)
 *
 * Problem solved
 * --------------
 *   -Δu + ∇p = f   in Ω
 *        div u = 0 in Ω
 *             u = u_exact on ∂Ω
 *
 * A Lagrange multiplier λ enforces the zero-mean pressure constraint:
 *
 *        ∫Ω p = 0
 *
 * The discrete system is
 *
 *   (∇u, ∇v) - (p, div v) + (div u, q) + (λ, q) + (p, μ) = (f, v)
 *
 * where μ is the test function associated with λ.
 *
 * Usage
 * -----
 * Compile the Rodin examples and run:
 *
 *   ./PETSc_Stokes
 *
 * PETSc solver options can be passed from the command line:
 *
 *   ./PETSc_Stokes -ksp_type gmres -pc_type gamg -ksp_rtol 1e-8 -ksp_monitor
 *
 * Output
 * ------
 * The program writes:
 *
 *   Stokes.mesh
 *   Stokes_velocity.gf
 *   Stokes_pressure.gf
 *   Stokes_velocity_exact.gf
 *
 * These can be visualized with the Rodin tools or external viewers
 * supporting the gridfunction format.
 */

#include <Rodin/PETSc.h>

#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

#include <petscksp.h>

#include <chrono>
#include <iostream>

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  // Mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, { 8, 8, 8 });
  mesh.scale(1.0 / (8 - 1));
  mesh.getConnectivity().compute(2, 3);
  mesh.getConnectivity().compute(3, 2);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);

  // FE spaces
  H1  uh(std::integral_constant<size_t, 2>{}, mesh, mesh.getSpaceDimension());
  H1  ph(std::integral_constant<size_t, 1>{}, mesh);
  P0g lh(mesh);

  std::cout << "Velocity FES size: " << uh.getSize() << "\n";
  std::cout << "Pressure FES size: " << ph.getSize() << "\n";
  std::cout << "LM       FES size: " << lh.getSize() << "\n";

  const auto pi = Constants::pi();

  // Manufactured solution
  VectorFunction u_exact{
    Sin(pi * F::x) * Cos(pi * F::y) * Cos(pi * F::z),
   -Cos(pi * F::x) * Sin(pi * F::y) * Cos(pi * F::z),
    Zero()
  };

  auto p_exact =
      Cos(2 * pi * F::x)
    * Cos(2 * pi * F::y)
    * Cos(2 * pi * F::z);

  VectorFunction f{
     3 * pi * pi * Sin(pi * F::x) * Cos(pi * F::y) * Cos(pi * F::z)
   - 2 * pi * Sin(2 * pi * F::x) * Cos(2 * pi * F::y) * Cos(2 * pi * F::z),

    -3 * pi * pi * Cos(pi * F::x) * Sin(pi * F::y) * Cos(pi * F::z)
   - 2 * pi * Cos(2 * pi * F::x) * Sin(2 * pi * F::y) * Cos(2 * pi * F::z),

    -2 * pi * Cos(2 * pi * F::x) * Cos(2 * pi * F::y) * Sin(2 * pi * F::z)
  };

  {
    // Unknowns and tests
    PETSc::Variational::TrialFunction u(uh); u.setName("u");
    PETSc::Variational::TrialFunction p(ph); p.setName("p");
    PETSc::Variational::TrialFunction l(lh); l.setName("lambda");

    PETSc::Variational::TestFunction v(uh);
    PETSc::Variational::TestFunction q(ph);
    PETSc::Variational::TestFunction m(lh);

    // Stokes + mean-pressure constraint
    Problem stokes(u, p, l, v, q, m);
    stokes =
        Integral(Jacobian(u), Jacobian(v))
      - Integral(p, Div(v))
      + Integral(Div(u), q)
      + Integral(l, q)
      + Integral(p, m)
      - Integral(f, v)
      + DirichletBC(u, u_exact);

    std::cout << "Assembling...\n";
    const auto t0 = std::chrono::high_resolution_clock::now();
    stokes.assemble();
    const auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Assembly time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << " ms\n";

    std::cout << "Solving...\n";
    Solver::KSP(stokes).solve();

    // Save solution
    mesh.save("Stokes.mesh", IO::FileFormat::MFEM);
    u.getSolution().save("Stokes_velocity.gf", IO::FileFormat::MFEM);
    p.getSolution().save("Stokes_pressure.gf", IO::FileFormat::MFEM);

    // Save interpolated exact velocity for visualization
    PETSc::Variational::GridFunction u_e(uh);
    u_e = u_exact;
    u_e.save("Stokes_velocity_exact.gf", IO::FileFormat::MFEM);
  }

  PetscFinalize();
  return 0;
}
