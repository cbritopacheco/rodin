#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <array>

#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;
using namespace Rodin::IO;
using namespace Rodin::External;

// Physical and algorithmic parameters
constexpr Real alpha       = 1.0;
constexpr Real beta        = 2.0;
constexpr Real hmin       = 0.1;
constexpr Real hmax       = 0.4;
constexpr Real targetOmega = 0.15;
constexpr Integer maxStepHom = 20;
constexpr Real mTarget     = 0.5;

int main()
{
  // Output setup
  const std::string name = "MinOneState";
  std::ofstream data("data/" + name + ".txt");
  data << std::setprecision(12);

  // Mesh generation: unit square with triangular elements
  MMG::Mesh mesh;
  size_t width = 16;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { width, width } );
  mesh.scale(1.0 / (width - 1.0));
  mesh.save("Omega.mesh");

  // MMG::Optimizer().setHMax(hmax)
  //                 .setHMin(hmin)
  //                 .setGradation(1.3)
  //                 .setHausdorff(0.01)
  //                 .optimize(mesh);

  mesh.getConnectivity().compute(1, 2);

  // Domain measure and cell count

  const Real volume = mesh.getVolume();
  const size_t cellCount = mesh.getCellCount();

  // Finite element spaces P1 for state, P0 for parameters
  P1 vh(mesh);
  P0 qh(mesh);

  auto u = TrialFunction(vh);
  auto v = TestFunction(vh);
  auto p = TrialFunction(vh);
  auto q = TestFunction(vh);

  // Fields: material coefficients, design variables, fluxes, sensitivities
  GridFunction a11(qh);
  GridFunction a12(qh);
  GridFunction a22(qh);
  GridFunction theta(qh);
  GridFunction sigma1(qh);
  GridFunction sigma2(qh);
  GridFunction tau1(qh);
  GridFunction tau2(qh);
  GridFunction M11(qh);
  GridFunction M12(qh);
  GridFunction M22(qh);
  GridFunction Case(qh);
  GridFunction compliance(qh);
  GridFunction f(qh);
  GridFunction objective(vh);

  f = 1.0;

  // Initial design: theta=1 everywhere
  theta = 1.0;
  a11 = theta * alpha + (1 - theta) * beta;
  a12 = 0.0;
  a22 = theta * alpha + (1 - theta) * beta;

  Real lagr = 0.0;
  std::vector<Real> input(6), R(5);

  // Define state and adjoint problems
  Problem state01(u, v);
  state01 = Integral(a11 * Dx(u), Dx(v))
          + Integral(a12 * Dx(u), Dy(v))
          + Integral(a12 * Dy(u), Dx(v))
          + Integral(a22 * Dy(u), Dy(v))
          - Integral(1.0 * v)
          + DirichletBC(u, Zero());


  Problem adjun01(p, v);
  adjun01 = Integral(a11 * Dx(p), Dx(v))
          + Integral(a12 * Dx(p), Dy(v))
          + Integral(a12 * Dy(p), Dx(v))
          + Integral(a22 * Dy(p), Dy(v))
          - Integral(1.0 * v)
          + DirichletBC(p, Zero());

  // Homogenization loop
  for (int stepHom = 0; stepHom < maxStepHom; ++stepHom)
  {
    state01.assemble();
    CG(state01).solve();

    adjun01.assemble();
    CG(adjun01).solve();

    p.getSolution().save("p.gf");

    // Compute fluxes
    sigma1 = a11 * Dx(u.getSolution()) + a12 * Dy(u.getSolution());
    sigma2 = a12 * Dx(u.getSolution()) + a22 * Dy(u.getSolution());
    tau1 = a11 * Dx(p.getSolution()) + a12 * Dy(p.getSolution());
    tau2 = a12 * Dx(p.getSolution()) + a22 * Dy(p.getSolution());

    // Sensitivity tensors
    M11 = -sigma1 * tau1;
    M12 = -(sigma1 * tau2 + sigma2 * tau1) * 0.5;
    M22 = -sigma2 * tau2;

    a11.save("a11.gf");
    a12.save("a12.gf");
    a22.save("a22.gf");
    Case.save("Case.gf");
    theta.save("theta.gf");
    u.getSolution().save("u.gf");
    p.getSolution().save("p.gf");

    compliance = f * u.getSolution();
    objective = f * u.getSolution() + lagr + theta;

    // Compute energies
    const Real energy = Integral(u.getSolution());
    const Real LagrangeHom = Integral(objective);

    data << stepHom << " " << lagr << " "
         << energy << " " << LagrangeHom << " "
         << volume << "\n";
    data.flush();

    // Lagrange multiplier update via secant
    double lagr1 = 0.0, lagr2 = 1e-1;
    double m1, m2, m;

    // Evaluate m1
    for (size_t i = 0; i < cellCount; ++i)
    {
      // input = { 1.0 / beta, 1.0 / alpha, lagr1, M11[i], M12[i], M22[i] };
      // theta[i] = updateThetaN(input);
    }
    m1 = Integral(theta).compute() - mTarget * volume;

    // Evaluate m2
    for (size_t i = 0; i < cellCount; ++i)
    {
      // input = {1.0/beta, 1.0/alpha, lagr2, M11[i], M12[i], M22[i]};
      // theta[i] = updateThetaN(input);
    }
    m2 = Integral(theta).compute() - mTarget * volume;

    if (m1 * m2 > 0)
    {
      std::cout << "Invalid bracket for lagrange: "
                << lagr1 << "(m1=" << m1 << "), "
                << lagr2 << "(m2=" << m2 << ")\n";
      break;
    }

    // Secant step
    lagr = (lagr1 * m2 - lagr2 * m1) / (m2 - m1);
    for (size_t i = 0; i < cellCount; ++i)
    {
      // input = {1.0/beta, 1.0/alpha, lagr, M11[i], M12[i], M22[i]};
      // theta[i] = updateThetaN(input);
    }
    m = Integral(theta).compute() - mTarget * volume;

    // Refine lagrange if needed
    for (int control = 0; control < 20; ++control)
    {
      if (std::abs(m) * 1e6 < hmax) break;
      if (m1 * m < 0)
      {
        lagr2 = lagr;
        m2 = m;
      }
      else
      {
        lagr1 = lagr; m1 = m;
      }
      lagr = (lagr1 * m2 - lagr2 * m1) / (m2 - m1);
      for (size_t i = 0; i < cellCount; ++i)
      {
        // input = {1.0/beta, 1.0/alpha, lagr, M11[i], M12[i], M22[i]};
        // theta[i] = updateThetaN(input);
      }
      m = Integral(theta).compute() - mTarget * volume;
    }

    // Update design and material fields
    for (size_t i = 0; i < cellCount; ++i)
    {
      // input = {1.0/beta, 1.0/alpha, lagr, M11[i], M12[i], M22[i]};
      // R     = updateConductivityN(input);
      // theta[i] = R[0];
      // a11[i]   = R[1];
      // a12[i]   = R[2];
      // a22[i]   = R[3];
      // Case[i]  = R[4];
    }
  }

  return 0;
}
