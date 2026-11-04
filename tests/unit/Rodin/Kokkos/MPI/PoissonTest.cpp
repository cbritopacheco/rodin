/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file PoissonTest.cpp
 * @brief Distributed device assembly of PETSc problems against the host.
 *
 * Rank 0 builds and partitions a grid, every rank reloads its shard, and the
 * same problem is solved with the distributed host assembly and with the
 * distributed device assembly. The owned coefficients of both solutions must
 * agree to the solver tolerance. The solutions then enter a second
 * problem as coefficients.
 */

#include "Common.h"

namespace Rodin::Tests::Unit
{
  /// @brief A Poisson problem with lifted Dirichlet values matches the host.
  void checkPoisson(Polytope::Type type)
  {
    using namespace Rodin::Variational;
    Context::MPI ctx(*g_env, *g_world);
    const HostMesh host = distributeFromRoot(ctx, type, 16);
    const DeviceMesh mesh(DeviceContext(ctx), host);
    P1<Real, HostMesh> hfes(host);
    P1<Real, DeviceMesh> dfes(mesh);
    EXPECT_EQ(hfes.getSize(), dfes.getSize());

    PETSc::Variational::TrialFunction hu(hfes);
    PETSc::Variational::TestFunction hv(hfes);
    PETSc::Variational::TrialFunction du(dfes);
    PETSc::Variational::TestFunction dv(dfes);

    RealFunction hf([](const Geometry::Point& p) { return 1 + p.x() * p.y(); });
    RealFunction hg([](const Geometry::Point& p) { return 2 * p.x() - p.y() + 1; });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });
    RealFunction dg(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 * p.x() - p.y() + 1;
    });

    Problem hp(hu, hv);
    hp = Integral(Grad(hu), Grad(hv)) + Integral(hu, hv) - Integral(hf, hv) +
      DirichletBC(hu, hg);
    Problem dp(du, dv);
    dp = Integral(Grad(du), Grad(dv)) + Integral(du, dv) - Integral(df, dv) +
      DirichletBC(du, dg);

    solve(hp);
    solve(dp);
    expectSameSystem(hp, dp);
    expectSameVector(hu.getSolution().getData(), du.getSolution().getData());

    // A second solve reuses the cached pattern and coordinates.
    solve(dp);
    expectSameVector(hu.getSolution().getData(), du.getSolution().getData());

    // The solutions enter a second problem as coefficients, through their
    // values and their gradients, ghosts included.
    const auto& hs = hu.getSolution();
    const auto& ds = du.getSolution();
    PETSc::Variational::TrialFunction hw(hfes);
    PETSc::Variational::TestFunction hz(hfes);
    PETSc::Variational::TrialFunction dw(dfes);
    PETSc::Variational::TestFunction dz(dfes);
    Problem hq(hw, hz);
    hq = Integral(Grad(hw), Grad(hz)) + Integral(hs * hw, hz) - Integral(hs, hz) -
      Integral(Dot(Grad(hs), Grad(hz))) + DirichletBC(hw, hg);
    Problem dq(dw, dz);
    dq = Integral(Grad(dw), Grad(dz)) + Integral(ds * dw, dz) - Integral(ds, dz) -
      Integral(Dot(Grad(ds), Grad(dz))) + DirichletBC(dw, dg);
    solve(hq);
    solve(dq);
    expectSameSystem(hq, dq);
    expectSameVector(hw.getSolution().getData(), dw.getSolution().getData());
  }

  /// @brief A stabilized P1/P1 Stokes problem matches the host.
  void checkStokes(Polytope::Type type)
  {
    using namespace Rodin::Variational;
    Context::MPI ctx(*g_env, *g_world);
    const HostMesh host = distributeFromRoot(ctx, type, 12);
    const DeviceMesh mesh(DeviceContext(ctx), host);
    P1 hvh(host, 2);
    P1 hqh(host);
    P1 dvh(mesh, 2);
    P1 dqh(mesh);

    PETSc::Variational::TrialFunction hu(hvh);
    PETSc::Variational::TestFunction hv(hvh);
    PETSc::Variational::TrialFunction hp(hqh);
    PETSc::Variational::TestFunction hq(hqh);
    PETSc::Variational::TrialFunction du(dvh);
    PETSc::Variational::TestFunction dv(dvh);
    PETSc::Variational::TrialFunction dp(dqh);
    PETSc::Variational::TestFunction dq(dqh);

    const Real eps = 1e-2;
    VectorFunction hf{[](const Geometry::Point& p) { return p.x() * p.y(); }, Real(1)};
    VectorFunction df(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return p.x() * p.y(); },
      Real(1));

    Problem hs(hu, hp, hv, hq);
    hs = Integral(Jacobian(hu), Jacobian(hv)) - Integral(hp, Div(hv)) +
      Integral(Div(hu), hq) + Integral(eps * hp, hq) - Integral(hf, hv) +
      DirichletBC(hu, VectorFunction{0, 0});
    Problem ds(du, dp, dv, dq);
    ds = Integral(Jacobian(du), Jacobian(dv)) - Integral(dp, Div(dv)) +
      Integral(Div(du), dq) + Integral(eps * dp, dq) - Integral(df, dv) +
      DirichletBC(du, VectorFunction{0, 0});

    solveDirect(hs);
    solveDirect(ds);
    expectSameSystem(hs, ds);
    expectSameVector(hu.getSolution().getData(), du.getSolution().getData());
    expectSameVector(hp.getSolution().getData(), dp.getSolution().getData());
  }

  /**
   * @brief A Neumann problem whose mean is held by a global multiplier
   * matches the host; the multiplier is owned by rank 0 and shared by all.
   */
  void checkMean(Polytope::Type type)
  {
    using namespace Rodin::Variational;
    Context::MPI ctx(*g_env, *g_world);
    const HostMesh host = distributeFromRoot(ctx, type, 12);
    const DeviceMesh mesh(DeviceContext(ctx), host);
    P1 hvh(host);
    P0g hgh(host);
    P1 dvh(mesh);
    P0g dgh(mesh);

    PETSc::Variational::TrialFunction hu(hvh);
    PETSc::Variational::TestFunction hv(hvh);
    PETSc::Variational::TrialFunction hl(hgh);
    PETSc::Variational::TestFunction hm(hgh);
    PETSc::Variational::TrialFunction du(dvh);
    PETSc::Variational::TestFunction dv(dvh);
    PETSc::Variational::TrialFunction dl(dgh);
    PETSc::Variational::TestFunction dm(dgh);

    // The small multiplier block gives its row the diagonal entry the LU
    // factorization of PETSc requires.
    const Real eps = 1e-2;
    RealFunction hf([](const Geometry::Point& p) { return 1 + p.x() * p.y(); });
    RealFunction df(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return 1 + p.x() * p.y(); });

    Problem hn(hu, hl, hv, hm);
    hn = Integral(Grad(hu), Grad(hv)) + Integral(hl, hv) + Integral(hu, hm) -
      Integral(eps * hl, hm) - Integral(hf, hv);
    Problem dn(du, dl, dv, dm);
    dn = Integral(Grad(du), Grad(dv)) + Integral(dl, dv) + Integral(du, dm) -
      Integral(eps * dl, dm) - Integral(df, dv);

    solveDirect(hn);
    solveDirect(dn);
    expectSameSystem(hn, dn);
    expectSameVector(hu.getSolution().getData(), du.getSolution().getData());
    expectSameVector(hl.getSolution().getData(), dl.getSolution().getData());
  }

  /**
   * @brief A quadratic Poisson problem matches the host, and so does a
   * second problem reading its solution, edge degrees of freedom included.
   */
  void checkQuadratic(Polytope::Type type)
  {
    using namespace Rodin::Variational;
    Context::MPI ctx(*g_env, *g_world);
    const HostMesh host = distributeFromRoot(ctx, type, 9);
    const DeviceMesh mesh(DeviceContext(ctx), host);
    const std::integral_constant<size_t, 2> k;
    H1 hfes(k, host);
    H1 dfes(k, mesh);
    EXPECT_EQ(hfes.getSize(), dfes.getSize());
    EXPECT_EQ(hfes.getShard().getSize(), dfes.getLocalSize());

    PETSc::Variational::TrialFunction hu(hfes);
    PETSc::Variational::TestFunction hv(hfes);
    PETSc::Variational::TrialFunction du(dfes);
    PETSc::Variational::TestFunction dv(dfes);

    RealFunction hf([](const Geometry::Point& p) { return 1 + p.x() * p.y(); });
    RealFunction hg([](const Geometry::Point& p) { return p.x() * p.x() - p.y(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });
    RealFunction dg(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return p.x() * p.x() - p.y();
    });

    Problem hp(hu, hv);
    hp = Integral(Grad(hu), Grad(hv)) + Integral(hu, hv) - Integral(hf, hv) +
      DirichletBC(hu, hg);
    Problem dp(du, dv);
    dp = Integral(Grad(du), Grad(dv)) + Integral(du, dv) - Integral(df, dv) +
      DirichletBC(du, dg);
    solve(hp);
    solve(dp);
    expectSameSystem(hp, dp);
    expectSameVector(hu.getSolution().getData(), du.getSolution().getData());

    const auto& hs = hu.getSolution();
    const auto& ds = du.getSolution();
    PETSc::Variational::TrialFunction hw(hfes);
    PETSc::Variational::TestFunction hz(hfes);
    PETSc::Variational::TrialFunction dw(dfes);
    PETSc::Variational::TestFunction dz(dfes);
    Problem hq(hw, hz);
    hq = Integral(Grad(hw), Grad(hz)) + Integral(hs * hw, hz) -
      Integral(Dot(Grad(hs), Grad(hz))) + DirichletBC(hw, hg);
    Problem dq(dw, dz);
    dq = Integral(Grad(dw), Grad(dz)) + Integral(ds * dw, dz) -
      Integral(Dot(Grad(ds), Grad(dz))) + DirichletBC(dw, dg);
    solve(hq);
    solve(dq);
    expectSameSystem(hq, dq);
    expectSameVector(hw.getSolution().getData(), dw.getSolution().getData());
  }

  /// @brief A regularized Taylor-Hood Stokes problem matches the host.
  void checkTaylorHood(Polytope::Type type)
  {
    using namespace Rodin::Variational;
    Context::MPI ctx(*g_env, *g_world);
    const HostMesh host = distributeFromRoot(ctx, type, 9);
    const DeviceMesh mesh(DeviceContext(ctx), host);
    const std::integral_constant<size_t, 2> k;
    H1 hvh(k, host, 2);
    P1 hqh(host);
    H1 dvh(k, mesh, 2);
    P1 dqh(mesh);
    EXPECT_EQ(hvh.getSize(), dvh.getSize());

    PETSc::Variational::TrialFunction hu(hvh);
    PETSc::Variational::TestFunction hv(hvh);
    PETSc::Variational::TrialFunction hp(hqh);
    PETSc::Variational::TestFunction hq(hqh);
    PETSc::Variational::TrialFunction du(dvh);
    PETSc::Variational::TestFunction dv(dvh);
    PETSc::Variational::TrialFunction dp(dqh);
    PETSc::Variational::TestFunction dq(dqh);

    const Real eps = 1e-2;
    VectorFunction hf{[](const Geometry::Point& p) { return p.x() * p.y(); }, Real(1)};
    VectorFunction df(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return p.x() * p.y(); },
      Real(1));

    Problem hs(hu, hp, hv, hq);
    hs = Integral(Jacobian(hu), Jacobian(hv)) - Integral(hp, Div(hv)) +
      Integral(Div(hu), hq) + Integral(eps * hp, hq) - Integral(hf, hv) +
      DirichletBC(hu, VectorFunction{0, 0});
    Problem ds(du, dp, dv, dq);
    ds = Integral(Jacobian(du), Jacobian(dv)) - Integral(dp, Div(dv)) +
      Integral(Div(du), dq) + Integral(eps * dp, dq) - Integral(df, dv) +
      DirichletBC(du, VectorFunction{0, 0});

    solveDirect(hs);
    solveDirect(ds);
    expectSameSystem(hs, ds);
    expectSameVector(hu.getSolution().getData(), du.getSolution().getData());
    expectSameVector(hp.getSolution().getData(), dp.getSolution().getData());
  }

  /**
   * @brief Distributed forms match the host, square or not, and enter
   * problems as preassembled forms, in a single block or in a mixed one.
   */
  void checkForms(Polytope::Type type)
  {
    using namespace Rodin::Variational;
    using Body = ProblemBody<::Mat, ::Vec, PetscScalar>;
    Context::MPI ctx(*g_env, *g_world);
    const HostMesh host = distributeFromRoot(ctx, type, 12);
    const DeviceMesh mesh(DeviceContext(ctx), host);
    P1 hfes(host);
    P1 dfes(mesh);
    P1 hvh(host, 2);
    P1 dvh(mesh, 2);

    PETSc::Variational::TrialFunction hu(hfes);
    PETSc::Variational::TestFunction hv(hfes);
    PETSc::Variational::TrialFunction du(dfes);
    PETSc::Variational::TestFunction dv(dfes);
    PETSc::Variational::TrialFunction hw(hvh);
    PETSc::Variational::TestFunction hz(hvh);
    PETSc::Variational::TrialFunction dw(dvh);
    PETSc::Variational::TestFunction dz(dvh);

    RealFunction hf([](const Geometry::Point& p) { return 1 + p.x() * p.y(); });
    RealFunction hg([](const Geometry::Point& p) { return 2 * p.x() - p.y() + 1; });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });
    RealFunction dg(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 * p.x() - p.y() + 1;
    });

    PETSc::Variational::BilinearForm hm(hu, hv);
    hm = Integral(hu, hv);
    hm.assemble();
    PETSc::Variational::BilinearForm dm(du, dv);
    dm = Integral(du, dv);
    dm.assemble();
    expectSameOperator(hm.getOperator(), dm.getOperator());

    PETSc::Variational::LinearForm hl(hv);
    hl = Integral(hf, hv);
    hl.assemble();
    PETSc::Variational::LinearForm dl(dv);
    dl = Integral(df, dv);
    dl.assemble();
    expectSameVector(hl.getVector(), dl.getVector());

    PETSc::Variational::BilinearForm hb(hw, hv);
    hb = Integral(Div(hw), hv);
    hb.assemble();
    PETSc::Variational::BilinearForm db(dw, dv);
    db = Integral(Div(dw), dv);
    db.assemble();
    expectSameOperator(hb.getOperator(), db.getOperator());

    const BilinearFormBase<::Mat>& hmass = hm;
    const BilinearFormBase<::Mat>& dmass = dm;
    const LinearFormBase<::Vec>& hload = hl;
    const LinearFormBase<::Vec>& dload = dl;
    Problem hp(hu, hv);
    hp = Body(Integral(Grad(hu), Grad(hv))) + hmass - hload + DirichletBC(hu, hg);
    Problem dp(du, dv);
    dp = Body(Integral(Grad(du), Grad(dv))) + dmass - dload + DirichletBC(du, dg);
    solve(hp);
    solve(dp);
    expectSameSystem(hp, dp);
    expectSameVector(hu.getSolution().getData(), du.getSolution().getData());

    // The penalty and the load of the pressure come second in the system.
    const Real eps = 1e-2;
    PETSc::Variational::BilinearForm he(hu, hv);
    he = Integral(eps * hu, hv);
    he.assemble();
    PETSc::Variational::BilinearForm de(du, dv);
    de = Integral(eps * du, dv);
    de.assemble();
    const BilinearFormBase<::Mat>& hpenalty = he;
    const BilinearFormBase<::Mat>& dpenalty = de;
    Problem hs(hw, hu, hz, hv);
    hs = Body(Integral(Jacobian(hw), Jacobian(hz))) - Integral(hu, Div(hz)) +
      Integral(Div(hw), hv) + hpenalty - hload +
      DirichletBC(hw, VectorFunction{0, 0});
    Problem ds(dw, du, dz, dv);
    ds = Body(Integral(Jacobian(dw), Jacobian(dz))) - Integral(du, Div(dz)) +
      Integral(Div(dw), dv) + dpenalty - dload +
      DirichletBC(dw, VectorFunction{0, 0});
    solveDirect(hs);
    solveDirect(ds);
    expectSameSystem(hs, ds);
    expectSameVector(hw.getSolution().getData(), dw.getSolution().getData());
    expectSameVector(hu.getSolution().getData(), du.getSolution().getData());
  }

  TEST(Kokkos_MPI_Forms, Triangles)
  {
    checkForms(Polytope::Type::Triangle);
  }

  TEST(Kokkos_MPI_Forms, Quadrilaterals)
  {
    checkForms(Polytope::Type::Quadrilateral);
  }

  TEST(Kokkos_MPI_Quadratic, Triangles)
  {
    checkQuadratic(Polytope::Type::Triangle);
  }

  TEST(Kokkos_MPI_Quadratic, Quadrilaterals)
  {
    checkQuadratic(Polytope::Type::Quadrilateral);
  }

  TEST(Kokkos_MPI_TaylorHood, Triangles)
  {
    checkTaylorHood(Polytope::Type::Triangle);
  }

  TEST(Kokkos_MPI_TaylorHood, Quadrilaterals)
  {
    checkTaylorHood(Polytope::Type::Quadrilateral);
  }

  TEST(Kokkos_MPI_Mean, Triangles)
  {
    checkMean(Polytope::Type::Triangle);
  }

  TEST(Kokkos_MPI_Mean, Quadrilaterals)
  {
    checkMean(Polytope::Type::Quadrilateral);
  }

  TEST(Kokkos_MPI_Stokes, Triangles)
  {
    checkStokes(Polytope::Type::Triangle);
  }

  TEST(Kokkos_MPI_Stokes, Quadrilaterals)
  {
    checkStokes(Polytope::Type::Quadrilateral);
  }

  TEST(Kokkos_MPI_Poisson, Triangles)
  {
    checkPoisson(Polytope::Type::Triangle);
  }

  TEST(Kokkos_MPI_Poisson, Quadrilaterals)
  {
    checkPoisson(Polytope::Type::Quadrilateral);
  }
}

/**
 * @brief Runs tests within the lifetimes of the MPI, Kokkos and PETSc runtimes.
 * @param argc Number of command-line arguments.
 * @param argv Command-line argument array.
 * @returns GoogleTest exit status.
 */
int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  Rodin::Tests::Unit::g_env = &env;
  Rodin::Tests::Unit::g_world = &world;
  ::testing::InitGoogleTest(&argc, argv);
  int res = 0;
  {
    ::Kokkos::ScopeGuard runtime(argc, argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    res = RUN_ALL_TESTS();
    PetscFinalize();
  }
  return res;
}
