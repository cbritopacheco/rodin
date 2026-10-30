/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [shape-source]
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/Solid.h>
#include <Rodin/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

using FES = P1<Math::SpatialVector<Real>, Mesh<Context::Local>>;

// Define boundary attributes
static constexpr Geometry::Attribute Gamma0 = 1; // Traction free boundary
static constexpr Geometry::Attribute GammaD = 2; // Homogenous Dirichlet boundary
static constexpr Geometry::Attribute GammaN = 3; // Inhomogenous Neumann boundary

// Lamé coefficients
static constexpr Real mu    = 0.3846;
static constexpr Real lambda = 0.5769;

// Optimization parameters
static constexpr size_t maxIt = 50;
static constexpr Real hmax  = 0.1;
static constexpr Real hmin  = 0.1 * hmax;
static constexpr Real ell  = 5.0;
static constexpr Real alpha = 4 * hmax;


template <class Data>
Real compliance(const GridFunction<FES, Data>& w)
{
  auto& vh = w.getFiniteElementSpace();
  TrialFunction u(vh);
  TestFunction  v(vh);
  BilinearForm  bf(u, v);
  bf = LinearElasticityIntegral(u, v)(lambda, mu);
  bf.assemble();
  return bf(w, w);
};

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <cantilever.mfem.mesh>" << std::endl;
    return EXIT_FAILURE;
  }
  const char* meshFile = argv[1];

  // Load mesh
  MMG::Mesh Omega;
  Omega.load(meshFile, IO::FileFormat::MFEM);
  Omega.save("Omega0.mesh", IO::FileFormat::MFEM);

  Alert::Info() << "Saved initial mesh to Omega0.mesh" << Alert::Raise;

  std::ofstream history("SimpleCantilever.csv");
  history.precision(17);
  history << "iteration,compliance,area,objective\n";
  IO::XDMF xdmf("SimpleCantilever");

  // Optimization loop
  for (size_t i = 0; i < maxIt; i++)
  {
    Omega.getConnectivity().compute(1, 2);
    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    // Finite element spaces
    int d = Omega.getSpaceDimension();
    P1 sh(Omega);
    P1 vh(Omega, d);

    // Pull-down force
    VectorFunction f{0, -1};

    //! [shape-state]
    // Elasticity equation
    TrialFunction u(vh);
    TestFunction  v(vh);
    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - BoundaryIntegral(f, v).over(GammaN)
               + DirichletBC(u, VectorFunction{0, 0}).on(GammaD);
    Solver::SparseLU stateSolver(elasticity);
    stateSolver.solve();
    if (!stateSolver.success())
      return EXIT_FAILURE;

    //! [shape-state]

    //! [shape-gradient]
    // Hilbert extension-regularization procedure
    TrialFunction g(vh);
    TestFunction  w(vh);

    auto e = 0.5 * (Jacobian(u.getSolution()) + Jacobian(u.getSolution()).T());
    auto Ae = 2.0 * mu * e + lambda * Trace(e) * IdentityMatrix(d);

    Problem hilbert(g, w);
    hilbert = Integral(alpha * Jacobian(g), Jacobian(w))
            + Integral(g, w)
            - BoundaryIntegral(Dot(Ae, e) - ell, Dot(BoundaryNormal(Omega), w)).over(Gamma0)
            + DirichletBC(g, VectorFunction{0, 0}).on({GammaD, GammaN});
    Solver::SparseLU gradientSolver(hilbert);
    gradientSolver.solve();
    if (!gradientSolver.success())
      return EXIT_FAILURE;

    //! [shape-gradient]

    g.getSolution().save("Theta.gf", IO::FileFormat::MFEM);
    Omega.save("Theta.mesh", IO::FileFormat::MFEM);

    const auto& dJ = g.getSolution();

    // Update objective
    const Real stateCompliance = compliance(u.getSolution());
    const Real area = Omega.getArea();
    const Real objective = stateCompliance + ell * area;
    if (!std::isfinite(objective))
      return EXIT_FAILURE;
    history << i << ',' << stateCompliance << ',' << area << ',' << objective << '\n';
    auto grid = xdmf.grid();
    grid.clear();
    grid.setMesh(Omega, IO::XDMF::MeshPolicy::Transient)
      .add("displacement", u.getSolution());
    xdmf.write(Real(i)).flush();
    Alert::Info() << "   | Objective: " << objective << Alert::Raise;

    //! [shape-update]
    // Make the displacement
    GridFunction norm(sh);
    norm = Frobenius(dJ);
    const Real maximum = norm.max();
    if (!std::isfinite(maximum))
      return EXIT_FAILURE;
    if (maximum == 0)
      break;
    g.getSolution() *= hmin / maximum;
    Omega.displace(g.getSolution());

    // Refine the mesh using MMG
    MMG::Optimizer().setHMax(hmax).setHMin(hmin).optimize(Omega);

    //! [shape-update]

    // Save mesh
    Omega.save("Omega.mesh", IO::FileFormat::MFEM);
  }

  Alert::Info() << "Saved final mesh to Omega.mesh" << Alert::Raise;

  history.flush();
  return history ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [shape-source]
