/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @example PETSc/ShapeOptimization/LevelSetCantilever2D.cpp
 * @brief Level-set based cantilever shape optimization in 2D using PETSc and MMG.
 *
 * This example solves a compliance minimization problem with perimeter penalization
 * by combining:
 * - linear elasticity for the state equation,
 * - a Hilbert-regularized shape gradient,
 * - level set advection,
 * - signed-distance reinitialization,
 * - remeshing with MMG.
 *
 * At each iteration, the code:
 * 1. optimizes the current mesh,
 * 2. trims the exterior subdomain,
 * 3. solves the elasticity state equation on the current shape,
 * 4. computes a regularized shape gradient,
 * 5. advects the level set function,
 * 6. reconstructs the domain with MMG,
 * 7. updates and records the objective value.
 *
 * The objective written to @c obj.txt is
 * @f[
 *   J(\Omega) = \text{compliance}(\Omega) + \ell |\Omega|,
 * @f]
 * where @f$\ell@f$ is the area penalization coefficient.
 *
 * @section input Input mesh
 * The example expects the initial mesh
 * @code
 * ../resources/examples/ShapeOptimization/LevelSetCantilever2D.mfem.mesh
 * @endcode
 *
 * @section usage Usage
 * Build the examples, then run for instance:
 *
 * @code
 * OMP_NUM_THREADS=8 ./examples/PETSc/ShapeOptimization/PETSc_LevelSetCantilever2D \
 *   -ksp_type cg \
 *   -pc_type gamg \
 *   -ksp_rtol 1e-8 \
 *   -ksp_max_it 1000 \
 *   -ksp_converged_reason \
 *   -ksp_monitor \
 *   -mat_block_size 2 \
 *   -pc_gamg_threshold 0.01
 * @endcode
 *
 * PETSc options are forwarded to the linear solves used for:
 * - the elasticity state equation,
 * - the Hilbert extension / regularization problem.
 *
 * A common alternative is to test without monitoring:
 *
 * @code
 * OMP_NUM_THREADS=8 ./examples/PETSc/ShapeOptimization/PETSc_LevelSetCantilever2D \
 *   -ksp_type cg -pc_type gamg -ksp_rtol 1e-8
 * @endcode
 *
 * @section output Output files
 * During the optimization, the program writes intermediate files such as:
 * - @c Omega0.mesh       : initial optimized mesh,
 * - @c Omega.mesh        : current trimmed domain,
 * - @c State.mesh        : mesh used for the state solve,
 * - @c State.gf          : state displacement,
 * - @c dJ.mesh           : mesh for the shape gradient,
 * - @c dJ.gf             : regularized shape gradient,
 * - @c Distance.mesh     : mesh for the signed-distance function,
 * - @c Distance.gf       : signed-distance function,
 * - @c Advect.mesh       : mesh after advection,
 * - @c Advect.gf         : advected level set,
 * - @c out/Omega.<k>.mesh: reconstructed domain at iteration @c k,
 * - @c obj.txt           : objective history.
 *
 * @section notes Notes
 * - The example uses PETSc for the linear algebra back-end.
 * - MMG is used both for mesh optimization and for level-set based remeshing.
 * - If remeshing fails at some iteration, the code reduces @c hmax and retries.
 * - The directory @c out/ should exist before running if you want to keep the
 *   saved iteration meshes.
 */

#include <Rodin/MMG.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Advection/Lagrangian.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

using FES = VectorP1<Mesh<Context::Local>>;

// Define interior and exterior for level set discretization
static constexpr Attribute interior = 1, exterior = 2;

// Define boundary attributes
static constexpr Attribute Gamma0 = 1, GammaD = 2, GammaN = 3, Gamma = 4;

// Lamé coefficients
static constexpr Real mu = 0.3846;
static constexpr Real lambda = 0.5769;

// Optimization parameters
static size_t maxIt = 300;
static Real hmax0 = 0.01;
static Real hmax = 0.01;
static Real hmin = 0.1 * hmax;
static Real hausd = 0.5 * hmin;
static Real ell = 0.4;
const Real dt = 0.5 * (hmax - hmin);
static Real alpha = 0.1;

// Compliance
template <class GridFunctionType>
Real compliance(const GridFunctionType& w)
{
  auto& vh = w.getFiniteElementSpace();
  PETSc::Variational::TrialFunction u(vh);
  PETSc::Variational::TestFunction  v(vh);
  BilinearForm  bf(u, v);
  bf = LinearElasticityIntegral(u, v)(lambda, mu);
  bf.assemble();
  return bf(w, w);
};

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  const char* meshFile = "../resources/examples/ShapeOptimization/LevelSetCantilever2D.mfem.mesh";

  // Load mesh
  MMG::Mesh th;
  th.load(meshFile, IO::FileFormat::MFEM);

  MMG::Optimizer().setHMax(hmax).setHMin(hmin).optimize(th);

  th.save("Omega0.mesh", IO::FileFormat::MEDIT);
  Alert::Info() << "Saved initial mesh to Omega0.mesh" << Alert::Raise;

  // Optimization loop
  std::vector<double> obj;
  std::ofstream fObj("obj.txt");
  size_t i = 0;
  while (i < maxIt)
  {
    try
    {
      MMG::Optimizer().setHMax(hmax)
                      .setHMin(hmin)
                      .setHausdorff(hausd)
                      .setAngleDetection(false)
                      .optimize(th);

      hmax = hmax0;
      hmin = 0.1 * hmax;
    }
    catch (const Alert::Exception& e)
    {
      hmax /= 2;
      hmin = 0.1 * hmax;
      Alert::Warning() << "Meshing failed at iteration " << i
        << ". Reducing hmax to " << hmax << " and retrying." << Alert::Raise;
      continue;
    }

    th.getConnectivity().compute(1, 2);
    th.getConnectivity().compute(0, 0);

    Alert::Info() << "----- Iteration: " << i << Alert::Raise;

    Alert::Info() << "   | Trimming mesh." << Alert::Raise;
    SubMesh trimmed = th.trim(exterior);
    trimmed.save("Omega.mesh", IO::FileFormat::MFEM);

    Alert::Info() << "   | Building finite element spaces." << Alert::Raise;
    const size_t d = th.getSpaceDimension();
    P1 sh(th);
    P1 vh(th, d);

    P1 shInt(trimmed);
    P1 vhInt(trimmed, d);

    Alert::Info() << "   | Solving state equation." << Alert::Raise;
    auto f = VectorFunction{0, -1};
    PETSc::Variational::TrialFunction u(vhInt);
    PETSc::Variational::TestFunction  v(vhInt);

    // Elasticity equation
    Problem elasticity(u, v);
    elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
               - BoundaryIntegral(f, v).over(GammaN)
               + DirichletBC(u, VectorFunction{0, 0}).on(GammaD);
    Solver::KSP(elasticity).solve();

    u.getSolution().save("State.gf", IO::FileFormat::MFEM);
    trimmed.save("State.mesh", IO::FileFormat::MFEM);

    Alert::Info() << "   | Computing shape gradient." << Alert::Raise;
    auto jac = Jacobian(u.getSolution());
    jac.traceOf(interior);
    auto e = 0.5 * (jac + jac.T());
    auto Ae = 2.0 * mu * e + lambda * Trace(e) * IdentityMatrix(d);
    auto n = FaceNormal(th);
    n.traceOf(interior);

    // Hilbert extension-regularization procedure
    PETSc::Variational::TrialFunction g(vh);
    PETSc::Variational::TestFunction w(vh);
    Problem hilbert(g, w);
    hilbert = Integral(alpha * alpha * Jacobian(g), Jacobian(w))
            + Integral(g, w)
            - FaceIntegral(Dot(Ae, e) - ell, Dot(n, w)).over(Gamma)
            + DirichletBC(g, VectorFunction{0, 0}).on(GammaN);
    Solver::KSP(hilbert).solve();

    auto& dJ = g.getSolution();
    dJ.save("dJ.gf", IO::FileFormat::MFEM);
    vh.getMesh().save("dJ.mesh", IO::FileFormat::MFEM);

    // Update objective
    double objective = compliance(u.getSolution()) + ell * th.getArea(interior);
    obj.push_back(objective);
    fObj << objective << "\n";
    fObj.flush();
    Alert::Info() << "   | Objective: " << obj.back() << Alert::Raise;
    Alert::Info() << "   | Distancing domain." << Alert::Raise;

    GridFunction dist(sh);
    Distance::Eikonal(dist).setInterior(interior)
                           .setInterface(Gamma)
                           .solve()
                           .sign();

    th.save("Distance.mesh", IO::FileFormat::MFEM);
    dist.save("Distance.gf", IO::FileFormat::MFEM);

    // Advect the level set function
    Alert::Info() << "   | Advecting the distance function." << Alert::Raise;
    GridFunction norm(sh);
    norm = Frobenius(dJ);
    dJ /= norm.max();

    TrialFunction advect(sh);
    TestFunction test(sh);

    Advection::Lagrangian(advect, test, dist, dJ).step(dt);

    th.save("Advect.mesh", IO::FileFormat::MFEM);
    advect.getSolution().save("Advect.gf", IO::FileFormat::MFEM);

    // Recover the implicit domain
    Alert::Info() << "   | Meshing the domain." << Alert::Raise;

    try
    {
      th = MMG::LevelSetDiscretizer().split(interior, {interior, exterior})
                                      .split(exterior, {interior, exterior})
                                      .setRMC(1e-6)
                                      .setHMax(hmax)
                                      .setHMin(hmin)
                                      .setHausdorff(hausd)
                                      .setAngleDetection(false)
                                      .setBoundaryReference(Gamma)
                                      .setBaseReferences(GammaD)
                                      .discretize(advect.getSolution());

      hmax = hmax0;
      hmin = 0.1 * hmax;
    }
    catch (const Alert::Exception& e)
    {
      hmax /= 2;
      hmin = 0.1 * hmax;
      Alert::Warning() << "Meshing failed at iteration " << i
        << ". Reducing hmax to " << hmax << " and retrying." << Alert::Raise;
      continue;
    }

    i++;
    th.save("out/Omega." + std::to_string(i) + ".mesh", IO::FileFormat::MEDIT);
  }

  Alert::Success() << "Saved final mesh to Omega.mesh" << Alert::Raise;

  PetscFinalize();

  return 0;
}

