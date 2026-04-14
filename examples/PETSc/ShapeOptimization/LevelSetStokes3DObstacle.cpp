/*
 * Single-file rewrite of the Stokes shape-optimization prototype.
 * Semantic choice is fixed here:
 *   - obstacle cells: attribute 2
 *   - fluid cells:    attribute 3
 *   - shape interface: attribute 13
 * The optimization target is the obstacle shape.
 *
 * OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
 * ./examples/PETSc/ShapeOptimization/LevelSetStokes3DObstacle \
 * -ksp_type preonly \
 * -pc_type lu \
 *  -ksp_converged_reason \
 *  -ksp_monitor \
 *  -pc_factor_shift_type nonzero \
 *  -pc_factor_shift_amount 1e-10 \
 *  -pc_factor_mat_solver_type mumps \
 *  -mat_mumps_icntl_7 7
 */

#include "Rodin/Geometry/Types.h"
#include <Rodin/MMG.h>
#include <Rodin/PETSc.h>
#include <Rodin/IO/XDMF.h>

#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Advection/Lagrangian.h>

#include <petscsys.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

using MeshType = Rodin::MMG::Mesh;

// Phase labels
static constexpr Attribute Obstacle   = 2;
static constexpr Attribute Fluid      = 3;
static constexpr Attribute GammaShape = 13;

// Stokes viscosity
static constexpr Real mu = 1.0;

// Optimization parameters
static constexpr size_t defaultMaxIt   = 60;
static constexpr Real   defaultHMax    = 0.2;
static constexpr Real   defaultAlpha   = 0.2;
static constexpr Real   defaultDt      = 0.5 * (defaultHMax - 0.1 * defaultHMax);
static constexpr Real   regularization = 1e-12;

static std::array<Real, 3> basis(int axis)
{
  if (axis < 0 || axis > 2)
    throw std::runtime_error("Axis must be 0, 1 or 2.");

  std::array<Real, 3> e{0.0, 0.0, 0.0};
  e[static_cast<size_t>(axis)] = 1.0;
  return e;
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  try
  {
    const char* meshFile = "sphere_init_aligned.mesh";
    const size_t maxIt   = defaultMaxIt;

    // Objective component K_ij, default K_11 in 0-based storage -> axis 1,1 as in original file.
    const int iAxis = 1;
    const int jAxis = 1;

    const Real hmax0 = defaultHMax;

    Real hmax  = hmax0;
    Real hmin  = 0.1 * hmax;
    Real hausd = 0.1 * hmin;

    const Real k = 0.1;
    const Real dt    = k * (defaultHMax - hmin);
    const Real alpha = 16 * (defaultHMax - hmin);

    // Load and pre-optimize initial mesh
    MeshType th;
    th.load(meshFile, IO::FileFormat::MEDIT);

    MMG::Optimizer()
      .setHMax(hmax)
      .setHMin(hmin)
      .setHausdorff(hausd)
      .optimize(th);

    th.save("Omega0.mesh", IO::FileFormat::MEDIT);

    // XDMF output
    IO::XDMF xdmf("out/LevelSetStokes3DObstacle");

    auto domainGrid = xdmf.grid("domain");
    domainGrid.setMesh(th, IO::XDMF::MeshPolicy::Transient);

    auto stateGrid = xdmf.grid("state");

    // Histories
    std::ofstream fObj("obj.txt");
    std::ofstream fVol("vol.txt");
    std::ofstream fAL("al.txt");

    if (!fObj || !fVol || !fAL)
      throw std::runtime_error("Failed to open output history files.");

    // Augmented Lagrangian parameters
    double lambdaAL = 0.0;
    double penalty  = 100;
    const double targetObstacleVolume = th.getVolume(Obstacle);

    // Objective basis vectors
    const auto ei = basis(iAxis);
    const auto ej = basis(jAxis);

    // Initial connectivity
    {
      auto& conn = th.getConnectivity();
      conn.compute(2, 3);
      conn.compute(3, 2);
      conn.compute(2, 1);
      conn.compute(1, 0);
      conn.compute(0, 0);
    }

    // Detect the fixed exterior boundary attributes once from the initial mesh
    FlatSet<Attribute> baseOuterBdr;
    for (auto it = th.getFace(); !it.end(); ++it)
    {
      const auto& face = *it;
      if (!face.isBoundary())
        continue;

      const auto a = face.getAttribute();
      if (!a || *a == GammaShape)
        continue;

      baseOuterBdr.insert(*a);
    }

    for (size_t it = 0; it < maxIt; ++it)
    {
      Alert::Info() << "----- Iteration: " << it << Alert::Raise;

      // ----------------------------------------------------------------------
      // Optimize current mesh
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Optimizing the domain." << Alert::Raise;
      try
      {
        MMG::Optimizer()
          .setHMax(hmax)
          .setHMin(hmin)
          .setHausdorff(hausd)
          .setAngleDetection(false)
          .optimize(th);

        hmax  = hmax0;
        hmin  = 0.1 * hmax;
        hausd = 0.1 * hmin;
      }
      catch (const Alert::Exception&)
      {
        hmax  /= 2;
        hmin   = 0.33 * hmax;
        hausd  = 0.1 * hmin;

        Alert::Warning()
          << "Mesh optimization failed at iteration " << it
          << ". Reducing hmax to " << hmax
          << " and retrying."
          << Alert::Raise;
        continue;
      }

      // ----------------------------------------------------------------------
      // Connectivity
      // ----------------------------------------------------------------------
      {
        auto& conn = th.getConnectivity();
        conn.compute(2, 3);
        conn.compute(3, 2);
        conn.compute(2, 1);
        conn.compute(1, 0);
        conn.compute(0, 0);
      }

      // ----------------------------------------------------------------------
      // Trim obstacle to get fluid mesh
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Trimming fluid mesh." << Alert::Raise;
      auto fluidMesh = th.trim(Obstacle);
      fluidMesh.save("Omega.mesh", IO::FileFormat::MEDIT);

      const size_t d = th.getSpaceDimension();

      // FE spaces
      auto velocitySpace = H1(std::integral_constant<size_t, 2>{}, fluidMesh, d);
      auto pressureSpace = H1(std::integral_constant<size_t, 1>{}, fluidMesh);
      P0g globalP0(fluidMesh);

      FlatSet<Attribute> fluidShapeBdr{GammaShape};
      FlatSet<Attribute> shapeInterface{GammaShape};

      // ----------------------------------------------------------------------
      // State Stokes solve
      // ----------------------------------------------------------------------
      PETSc::Variational::TrialFunction up(velocitySpace); up.setName("u");
      PETSc::Variational::TrialFunction pp(pressureSpace); pp.setName("p");

      PETSc::Variational::TestFunction vp(velocitySpace);
      PETSc::Variational::TestFunction qp(pressureSpace);

      Problem state(up, pp, vp, qp);
      state =
          Integral(mu * Jacobian(up), Jacobian(vp))
        - Integral(pp, Div(vp))
        + Integral(Div(up), qp)
        + regularization * Integral(pp, qp)
        + DirichletBC(up, VectorFunction{ej[0], ej[1], ej[2]}).on(fluidShapeBdr)
        + DirichletBC(up, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      Alert::Info() << "   | Assembling state Stokes." << Alert::Raise;
      state.assemble();
      state.setFieldSplits();

      Alert::Info() << "   | Solving state Stokes." << Alert::Raise;
      Solver::KSP(state)
        .setPrefix("state_")
        .solve();

      auto u = GridFunction(velocitySpace);
      auto p = GridFunction(pressureSpace);
      u = up.getSolution();
      p = pp.getSolution();

      // ----------------------------------------------------------------------
      // Adjoint Stokes solve
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Solving adjoint Stokes." << Alert::Raise;

      PETSc::Variational::TrialFunction ua(velocitySpace); ua.setName("v");
      PETSc::Variational::TrialFunction pa(pressureSpace); pa.setName("q");

      PETSc::Variational::TestFunction va(velocitySpace);
      PETSc::Variational::TestFunction qa(pressureSpace);

      Problem adj(ua, pa, va, qa);
      adj =
          Integral(mu * Jacobian(ua), Jacobian(va))
        - Integral(pa, Div(va))
        + Integral(Div(ua), qa)
        + regularization * Integral(pa, qa)
        + DirichletBC(ua, VectorFunction{ei[0], ei[1], ei[2]}).on(fluidShapeBdr)
        + DirichletBC(ua, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      Alert::Info() << "   | Assembling adjoint Stokes." << Alert::Raise;
      adj.assemble();
      adj.setFieldSplits();

      Alert::Info() << "   | Solving adjoint Stokes." << Alert::Raise;
      Solver::KSP(adj)
        .setPrefix("adj_")
        .solve();

      auto v = GridFunction(velocitySpace);
      auto q = GridFunction(pressureSpace);
      v = ua.getSolution();
      q = pa.getSolution();

      // ----------------------------------------------------------------------
      // Objective evaluation
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Evaluating objective." << Alert::Raise;

      auto ju = Jacobian(u);

      auto nObj = FaceNormal(th);
      nObj.traceOf(Obstacle);

      const auto eiVec = VectorFunction{ei[0], ei[1], ei[2]};
      auto sigma = mu * (ju + ju.T()) - p * IdentityMatrix(d);
      auto traction = sigma * nObj;

      P0g p0Obj(th);
      TestFunction z0Obj(p0Obj);
      LinearForm lfObj(z0Obj);
      lfObj = FaceIntegral(Dot(traction, eiVec), z0Obj).over(shapeInterface);
      lfObj.assemble();

      GridFunction oneObj(p0Obj);
      oneObj = 1.0;

      const double J = -lfObj(oneObj);

      Alert::Info() << "   | Objective: " << J << Alert::Raise;

      // ----------------------------------------------------------------------
      // Shape gradient
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Computing shape gradient." << Alert::Raise;

      auto juShape = Jacobian(u);
      auto jvShape = Jacobian(v);

      auto eu = 0.5 * (juShape + juShape.T());
      auto ev = 0.5 * (jvShape + jvShape.T());

      // Boundary integrand used for the Hilbert extension forcing
      auto G = 2.0 * mu * Dot(eu, ev);

      const double obstacleVolume = th.getVolume(Obstacle);
      const double violation = obstacleVolume - targetObstacleVolume;
      const double shift = -lambdaAL + penalty * violation;

      P0g p0Stats(th);
      TestFunction z0Stats(p0Stats);

      LinearForm lfArea(z0Stats);
      lfArea = FaceIntegral(z0Stats).over(shapeInterface);
      lfArea.assemble();

      LinearForm lfG(z0Stats);
      lfG = FaceIntegral(G, z0Stats).over(shapeInterface);
      lfG.assemble();

      GridFunction oneStats(p0Stats);
      oneStats = 1.0;

      const double shapeArea = lfArea(oneStats);
      const double meanG = (shapeArea > 0.0) ? lfG(oneStats) / shapeArea : 0.0;

      auto n = FaceNormal(th);
      n.traceOf(Obstacle);

      // ----------------------------------------------------------------------
      // Hilbert extension
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Solving vector Hilbert extension." << Alert::Raise;

      P1 vh(th, d);

      PETSc::Variational::TrialFunction theta(vh); theta.setName("theta");
      PETSc::Variational::TestFunction psi(vh);

      Problem hilbert(theta, psi);
      hilbert =
          Integral(alpha * alpha * Jacobian(theta), Jacobian(psi))
        + Integral(theta, psi)
        + FaceIntegral((G + shift) * Dot(n, psi)).over(shapeInterface)
        + DirichletBC(theta, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      Alert::Info() << "   | Assembling Hilbert extension." << Alert::Raise;
      hilbert.assemble();

      Alert::Info() << "   | Solving Hilbert extension." << Alert::Raise;
      Solver::KSP(hilbert)
        .setPrefix("hilbert_")
        .solve();

      // ----------------------------------------------------------------------
      // Signed-distance reconstruction
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Distancing obstacle phase." << Alert::Raise;

      P1 sh(th);
      GridFunction dist(sh);
      Distance::Eikonal(dist)
        .setInterior(Fluid)
        .setInterface(shapeInterface)
        .solve()
        .sign();

      // ----------------------------------------------------------------------
      // Advection
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Advecting level set." << Alert::Raise;

      auto& V = theta.getSolution();

      P1 shNorm(th);
      GridFunction speed(shNorm);
      speed = Frobenius(V);

      const Real vmax = speed.max();
      if (vmax > 1.0)
        V /= vmax;

      TrialFunction advect(sh);
      TestFunction test(sh);
      Advection::Lagrangian(advect, test, dist, V).step(dt);

      // ----------------------------------------------------------------------
      // XDMF snapshot before remeshing
      // ----------------------------------------------------------------------
      domainGrid.clear();
      domainGrid.setMesh(th, IO::XDMF::MeshPolicy::Transient);
      domainGrid.add("dist",   dist,                 IO::XDMF::Center::Node);
      domainGrid.add("theta",  theta.getSolution(),  IO::XDMF::Center::Node);
      domainGrid.add("speed",  speed,                IO::XDMF::Center::Node);
      domainGrid.add("advect", advect.getSolution(), IO::XDMF::Center::Node);

      stateGrid.clear();
      stateGrid.setMesh(fluidMesh, IO::XDMF::MeshPolicy::Transient);
      stateGrid.add("u", u, IO::XDMF::Center::Node);
      stateGrid.add("p", p, IO::XDMF::Center::Node);
      stateGrid.add("v", v, IO::XDMF::Center::Node);
      stateGrid.add("q", q, IO::XDMF::Center::Node);

      xdmf.write(it).flush();

      // ----------------------------------------------------------------------
      // Remeshing
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Remeshing." << Alert::Raise;

      try
      {
        th =
          MMG::LevelSetDiscretizer()
            .split(Fluid,    {Fluid, Obstacle})
            .split(Obstacle, {Fluid, Obstacle})
            .setRMC(1e-5)
            .setHMax(hmax)
            .setHMin(hmin)
            .setHausdorff(hausd)
            .setAngleDetection(false)
            .setBoundaryReference(GammaShape)
            .setBaseReferences(baseOuterBdr)
            .discretize(advect.getSolution());

        // NOTE: MMG adds a layer of old interface triangles with attribute
        // 2... We need to remove them for the optimization to work correctly.
        // This is a quirk of the current implementation and may be fixed in
        // future MMG versions.
        for (auto it = th.getFace(); it; ++it)
        {
          if (it->getAttribute() == Geometry::Attribute(Obstacle))
            th.setAttribute(it.key(), {});
        }
      }
      catch (const Alert::Exception&)
      {
        hmax  /= 2;
        hmin   = 0.33 * hmax;
        hausd  = 0.1 * hmin;

        Alert::Warning()
          << "Remeshing failed. Retrying with hmax=" << hmax
          << Alert::Raise;
        continue;
      }

      // ----------------------------------------------------------------------
      // Augmented Lagrangian update
      // ----------------------------------------------------------------------
      const double Jaug = J - lambdaAL * violation + 0.5 * penalty * violation * violation;

      fObj << J << " " << Jaug << "\n";
      fVol << obstacleVolume << " " << violation << "\n";
      fAL  << lambdaAL << " " << penalty << "\n";

      fObj.flush();
      fVol.flush();
      fAL.flush();

      lambdaAL += -penalty * violation;
      if (std::abs(violation) > 0.01 * targetObstacleVolume)
        penalty = std::min(1.0e6, 1.2 * penalty);

      Alert::Info()
        << "   | J=" << J
        << ", Jaug=" << Jaug
        << ", Vobs=" << obstacleVolume
        << ", c=" << violation
        << ", mean(G)=" << meanG
        << ", shift=" << shift
        << ", lambda=" << lambdaAL
        << ", penalty=" << penalty
        << Alert::Raise;

      th.save("out/Omega." + std::to_string(it) + ".mesh", IO::FileFormat::MEDIT);
    }

    xdmf.close();

    Alert::Success()
      << "Saved final mesh sequence to out/Omega.*.mesh"
      << Alert::Raise;

    PetscFinalize();
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << "LevelSetStokes3DObstacle failed: " << e.what() << '\n';
    PetscFinalize();
    return 1;
  }
}
