/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @example PETSc/ShapeOptimization/LevelSetCantilever3D.cpp
 * @brief 3D cantilever shape optimization using level sets, PETSc, and MMG.
 *
 * This example solves a compliance minimization problem with volume
 * penalization using a level-set representation of the domain.
 *
 * The algorithm combines:
 *  - linear elasticity for the state equation,
 *  - a Hilbert-regularized shape gradient,
 *  - signed-distance reconstruction of the level set,
 *  - semi-Lagrangian level set advection,
 *  - remeshing with MMG.
 *
 * At each optimization iteration the program:
 *  1. optimizes the current mesh with MMG,
 *  2. trims the exterior subdomain,
 *  3. solves the elasticity state equation,
 *  4. computes a regularized shape gradient,
 *  5. advects the level set function,
 *  6. reconstructs the domain using MMG,
 *  7. evaluates the objective.
 *
 * The objective recorded in @c obj.txt is
 *
 * @f[
 *   J(\Omega) = \text{compliance}(\Omega) + \ell |\Omega|
 * @f]
 *
 * where @f$\ell@f$ is the volume penalization parameter.
 *
 * @section input Input mesh
 *
 * The program expects the initial mesh:
 *
 * @code
 * ../resources/examples/ShapeOptimization/LevelSetCantilever3D.medit.mesh
 * @endcode
 *
 * This mesh contains:
 *  - interior and exterior region labels,
 *  - boundary attributes for Dirichlet and Neumann boundaries.
 *
 * @section usage Usage
 *
 * Build the Rodin examples and run for instance:
 *
 * @code
 * OMP_NUM_THREADS=8 ./examples/PETSc/ShapeOptimization/PETSc_LevelSetCantilever3D \
 *   -ksp_type cg \
 *   -pc_type gamg \
 *   -ksp_rtol 1e-8 \
 *   -ksp_max_it 1000 \
 *   -ksp_converged_reason \
 *   -ksp_monitor \
 *   -pc_gamg_threshold 0.01
 * @endcode
 *
 * PETSc options control the linear solvers used for:
 *  - the elasticity state equation,
 *  - the Hilbert regularization problem.
 *
 * @section output Output files
 *
 * During the optimization the following files are produced:
 *
 * - @c Omega0.mesh      : initial mesh
 * - @c Optimized.mesh   : mesh after MMG optimization
 * - @c Trimmed.mesh     : mesh restricted to the current domain
 * - @c State.mesh       : mesh used for the elasticity solve
 * - @c State.gf         : displacement field
 * - @c dJ.mesh          : mesh used for the shape gradient
 * - @c dJ.gf            : regularized shape gradient
 * - @c Distance.mesh    : mesh for the signed-distance field
 * - @c Distance.gf      : signed-distance function
 * - @c Advect.mesh      : mesh after advection
 * - @c Advect.gf        : advected level set
 * - @c Omega.mesh       : reconstructed domain
 * - @c obj.txt          : objective history
 *
 * @section notes Notes
 *
 * - PETSc provides the linear algebra backend.
 * - MMG is used for mesh optimization and level-set discretization.
 * - If remeshing fails at some iteration, the algorithm reduces @c hmax
 *   and retries with a finer mesh.
 * - For visualization, the mesh and grid functions can be opened with
 *   tools supporting the MEDIT format.
 */

#include <Rodin/MMG.h>
#include <Rodin/PETSc.h>

#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Advection/Lagrangian.h>

#include <fstream>
#include <vector>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// Define interior and exterior for level set discretization
static constexpr Attribute Interior = 3, Exterior = 2;

// Define boundary attributes
static constexpr Attribute GammaD = 4, GammaN = 7, Gamma = 10, Gamma0 = 3;

// Lamé coefficients
static constexpr Real mu = 0.3846;
static constexpr Real lambda = 0.5769;

// Optimization parameters
static constexpr size_t maxIt = 300;
static constexpr Real eps = 1e-12;
static constexpr Real hgrad = 1.6;
static constexpr Real ell = 0.1;
static Real elementStep = 0.5;
static Real hmax0 = 0.1;
static Real hmax = hmax0;
static Real hmin = 0.1 * hmax;
static Real hausd = 0.1 * hmin;
static size_t hmaxIt = maxIt / 2;
const Real k = 0.5;
const Real dt = k * (hmax - hmin) / 2;
static Real alpha = dt;

using FES = VectorP1<Mesh<Context::Local>>;

template <class GridFunctionType>
Real compliance(const GridFunctionType& w)
{
  auto& vh = w.getFiniteElementSpace();
  PETSc::Variational::TrialFunction u(vh);
  PETSc::Variational::TestFunction  v(vh);
  BilinearForm bf(u, v);
  bf = LinearElasticityIntegral(u, v)(lambda, mu);
  bf.assemble();
  return bf(w, w);
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  const char* meshFile = "../resources/examples/ShapeOptimization/LevelSetCantilever3D.medit.mesh";

  // Load mesh
  MMG::Mesh th;
  th.load(meshFile, IO::FileFormat::MEDIT);

  P1 sh(th);

  Alert::Info() << "Saved initial mesh to Omega0.mesh" << Alert::Raise;
  th.save("Omega0.mesh", IO::FileFormat::MEDIT);

  // Optimization loop
  std::vector<double> obj;
  std::ofstream fObj("obj.txt");
  size_t i = 0;
  while (i < maxIt)
  {
    Alert::Info() << "----- Iteration: " << i << Alert::Raise;
    Alert::Info() << "   | Optimizing the domain..." << Alert::Raise;

    try
    {
      MMG::Optimizer().setHMax(hmax)
                      .setHMin(hmin)
                      // .setGradation(1.1)
                      // .setHausdorff(hausd)
                      .setAngleDetection(false)
                      .optimize(th);

      hmax = hmax0;
      hmin = 0.1 * hmax;
    }
    catch (const Alert::Exception& e)
    {
      hmax /= 2;
      hmin = 0.1 * hmax;
      Alert::Warning() << "Mesh optimization failed at iteration " << i
                       << ". Reducing hmax to " << hmax
                       << " and retrying." << Alert::Raise;
      continue;
    }

    th.save("Optimized.mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "   | Computing connectivity." << Alert::Raise;
    auto& conn = th.getConnectivity();
    conn.discover(3, 2);
    conn.discover(3, 1);
    conn.restrict(1, 0);
    conn.restrict(2, 0);
    conn.restrict(2, 3);
    conn.discover(0, 0);

    Alert::Info() << "   | Trimming mesh." << Alert::Raise;
    SubMesh trimmed = th.trim(Exterior);
    trimmed.save("Trimmed.mesh", IO::FileFormat::MEDIT);

    Alert::Info() << "   | Building finite element spaces." << Alert::Raise;
    const size_t d = th.getSpaceDimension();

    P1 sh(th);
    P1 vh(th, d);

    P1 shInt(trimmed);
    P1 vhInt(trimmed, d);

    Alert::Info() << "   | Distancing domain." << Alert::Raise;
    GridFunction dist(sh);
    Distance::Eikonal(dist).setInterior(Interior)
                           .setInterface(Gamma)
                           .solve()
                           .sign();

    dist.getFiniteElementSpace().getMesh().save("Distance.mesh");
    dist.save("Distance.gf");

    Alert::Info() << "   | Solving state equation." << Alert::Raise;
    auto f = VectorFunction{0, 0, -1};

    PETSc::Variational::TrialFunction u(vhInt);
    PETSc::Variational::TestFunction  v(vhInt);

    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - BoundaryIntegral(f, v).over(GammaN)
               + DirichletBC(u, VectorFunction{0, 0, 0}).on(GammaD);
    Solver::KSP(elasticity).solve();

    u.getSolution().save("State.gf");
    trimmed.save("State.mesh");

    Alert::Info() << "   | Computing shape gradient." << Alert::Raise;
    auto jac = Jacobian(u.getSolution());
    jac.traceOf(Interior);

    auto e = 0.5 * (jac + jac.T());
    auto Ae = 2.0 * mu * e + lambda * Trace(e) * IdentityMatrix(d);

    auto n = FaceNormal(th);
    n.traceOf(Interior);

    PETSc::Variational::TrialFunction g(vh);
    PETSc::Variational::TestFunction  w(vh);

    Problem hilbert(g, w);
    hilbert = Integral(alpha * Jacobian(g), Jacobian(w))
            + Integral(g, w)
            - FaceIntegral(Dot(Ae, e) - ell, Dot(n, w)).over(Gamma)
            + DirichletBC(g, VectorFunction{0, 0, 0}).on(GammaN);
    Solver::KSP(hilbert).solve();

    auto& dJ = g.getSolution();
    vh.getMesh().save("dJ.mesh");
    dJ.save("dJ.gf");

    // Update objective
    double objective = compliance(u.getSolution()) + ell * th.getVolume(Interior);
    obj.push_back(objective);
    fObj << objective << "\n";
    fObj.flush();

    Alert::Info() << "   | Objective: " << obj.back() << Alert::Raise;

    // Advect the level set function
    Alert::Info() << "   | Advecting the distance function." << Alert::Raise;
    GridFunction norm(sh);
    norm = Frobenius(dJ);
    dJ /= norm.max();

    TrialFunction advect(sh);
    TestFunction  test(sh);

    Advection::Lagrangian(advect, test, dist, dJ).step(dt);

    th.save("Advect.mesh");
    advect.getSolution().save("Advect.gf");

    // Recover the implicit domain
    Alert::Info() << "   | Meshing the domain." << Alert::Raise;
    try
    {
      th = MMG::LevelSetDiscretizer().setHMax(hmax)
                                     .setHMin(hmin)
                                     // .setHausdorff(hausd)
                                     .setAngleDetection(false)
                                     .setRMC(1e-5)
                                     .setBaseReferences(GammaD)
                                     .setBoundaryReference(Gamma)
                                     .discretize(advect.getSolution());

      hmax = hmax0;
      hmin = 0.1 * hmax;
    }
    catch (const Alert::Exception& e)
    {
      hmax /= 2;
      hmin = 0.1 * hmax;
      Alert::Warning() << "Meshing failed at iteration " << i
                       << ". Reducing hmax to " << hmax
                       << " and retrying." << Alert::Raise;
      continue;
    }

    trimmed.save("out/Omega." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);

    i++;
  }

  Alert::Success() << "Saved final mesh to Omega.mesh" << Alert::Raise;

  PetscFinalize();

  return 0;
}
