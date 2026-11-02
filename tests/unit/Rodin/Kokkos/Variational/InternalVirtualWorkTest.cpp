/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solid.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>

#include "Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;

  constexpr std::integral_constant<size_t, 2> Quadratic{};

  constexpr std::integral_constant<size_t, 3> Cubic{};

  /**
   * @brief Sets a small displacement, keeping the deformation invertible.
   * @param data Coefficients.
   */
  template <class Data>
  void fillDisplacement(Data& data)
  {
    for (Index i = 0; i < static_cast<Index>(data.size()); ++i)
      data.coeffRef(i) = Real(0.005) * std::sin(Real(0.7) * i);
  }

  /// @brief The Neo-Hookean tangent and residual match the host.
  template <class HostFES, class DeviceFES>
  void checkForms(const HostFES& hfes, const DeviceFES& dfes)
  {
    GridFunction hw(hfes);
    GridFunction dw(dfes);
    fillDisplacement(hw.getData());
    fillDisplacement(dw.getData());
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);
    const Solid::NeoHookean law(Real(1.2), Real(0.8));
    const auto hivw = Solid::InternalVirtualWork(law, hw);
    const auto divw = Solid::InternalVirtualWork(law, dw);
    expectSameForm(hu, hv, hivw.Tangent(hu, hv), du, dv, divw.Tangent(du, dv));
    expectSameLinearForm(hv, hivw.Residual(hv), dv, divw.Residual(dv));
  }

  /// @brief Forms on vector P1 and quadratic spaces match the host.
  void checkForms(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    {
      SCOPED_TRACE("P1");
      checkForms(P1(source, d), P1(mesh, d));
    }
    {
      SCOPED_TRACE("H1 quadratic");
      checkForms(H1(Quadratic, source, d), H1(Quadratic, mesh, d));
    }
  }

  /// @brief A Newton solve under a body force converges to the host solution.
  void checkNewton(const HostMesh& source)
  {
    ASSERT_EQ(source.getSpaceDimension(), 2);
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source, 2);
    P1 dfes(mesh, 2);
    GridFunction hw(hfes);
    GridFunction dw(dfes);
    hw.getData().setZero();
    dw.getData().setZero();
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);
    const Solid::NeoHookean law(Real(1.2), Real(0.8));
    const auto hivw = Solid::InternalVirtualWork(law, hw);
    const auto divw = Solid::InternalVirtualWork(law, dw);

    Problem hp(hu, hv);
    hp = hivw(hu, hv) - Integral(VectorFunction{0, Real(-0.1)}, hv) +
      DirichletBC(hu, VectorFunction{0, 0}).on(1);
    Problem dp(du, dv);
    dp = divw(du, dv) - Integral(VectorFunction{0, Real(-0.1)}, dv) +
      DirichletBC(du, VectorFunction{0, 0}).on(1);

    Solver::SparseLU hs(hp);
    Solver::NewtonSolver hn(hs);
    hn.setMaxIterations(20).setAbsoluteTolerance(1e-12).setRelativeTolerance(1e-10);
    hn.solve(hw);
    Solver::SparseLU ds(dp);
    Solver::NewtonSolver dn(ds);
    dn.setMaxIterations(20).setAbsoluteTolerance(1e-12).setRelativeTolerance(1e-10);
    dn.solve(dw);

    ASSERT_TRUE(hn.converged());
    ASSERT_TRUE(dn.converged());
    EXPECT_EQ(hn.getReport().iterations, dn.getReport().iterations);
    ASSERT_GT(hw.getData().norm(), Real(0));
    EXPECT_LE((hw.getData() - dw.getData()).norm(), Real(1e-10) * hw.getData().norm());
  }

  TEST(Kokkos_Variational_InternalVirtualWork, FormsOnTriangles)
  {
    checkForms(withClosure(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_InternalVirtualWork, FormsOnTetrahedra)
  {
    checkForms(withClosure(makeTetrahedronMesh()));
  }

  TEST(Kokkos_Variational_InternalVirtualWork, FormsOnMixedMesh)
  {
    checkForms(withClosure(makeMixedMesh()));
  }

  /// One team per cell holds the local matrices beyond the thread capacity.
  TEST(Kokkos_Variational_InternalVirtualWork, FormsOfDegreeThreeOnHexahedra)
  {
    const auto source = withClosure(makeHexahedronMesh());
    const DeviceMesh mesh(TestContext{}, source);
    checkForms(H1(Cubic, source, 3), H1(Cubic, mesh, 3));
  }

  TEST(Kokkos_Variational_InternalVirtualWork, NewtonOnTriangles)
  {
    checkNewton(withClosure(makeTriangleMesh()));
  }
}
