#include <gtest/gtest.h>

#include "Rodin/Variational.h"
#include "Rodin/Assembly/Default.h"
#include "Rodin/Geometry/Mesh.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /// @brief Verifies mixed order grad grad assembles for variational H1 quadrature rule by checking exact expected values, form assembly.
  TEST(Rodin_Variational_H1QuadratureRule, MixedOrder_GradGrad_Assembles)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 0);
    H1<2, Real> fesTr(std::integral_constant<size_t, 2>{}, mesh);
    H1<1, Real> fesTe(std::integral_constant<size_t, 1>{}, mesh);

    TrialFunction u(fesTr);
    TestFunction v(fesTe);

    BilinearForm bf(u, v);
    bf = Integral(Grad(u), Grad(v));
    bf.assemble();

    const auto& mat = bf.getOperator();
    EXPECT_EQ(mat.rows(), static_cast<Eigen::Index>(fesTe.getSize()));
    EXPECT_EQ(mat.cols(), static_cast<Eigen::Index>(fesTr.getSize()));
    EXPECT_GT(mat.norm(), 0.0);
  }

  /// @brief Verifies mixed order jacobian assembles for variational H1 quadrature rule by checking exact expected values, form assembly.
  TEST(Rodin_Variational_H1QuadratureRule, MixedOrder_Jacobian_Assembles)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 0);
    const size_t vdim = 2;
    H1<2, Math::SpatialVector<Real>> fesTr(std::integral_constant<size_t, 2>{}, mesh, vdim);
    H1<1, Math::SpatialVector<Real>> fesTe(std::integral_constant<size_t, 1>{}, mesh, vdim);

    TrialFunction u(fesTr);
    TestFunction v(fesTe);

    BilinearForm bf(u, v);
    bf = Integral(Jacobian(u), Jacobian(v));
    bf.assemble();

    const auto& mat = bf.getOperator();
    EXPECT_EQ(mat.rows(), static_cast<Eigen::Index>(fesTe.getSize()));
    EXPECT_EQ(mat.cols(), static_cast<Eigen::Index>(fesTr.getSize()));
    EXPECT_GT(mat.norm(), 0.0);
  }

  /// @brief Verifies vector mass matches grid function linear form for variational H1 quadrature rule by checking tolerance-based numerical results, exact expected values, form assembly.
  TEST(Rodin_Variational_H1QuadratureRule, VectorMass_MatchesGridFunctionLinearForm)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 0);

    const size_t vdim = 2;
    H1<2, Math::SpatialVector<Real>> fes(std::integral_constant<size_t, 2>{}, mesh, vdim);
    TrialFunction u(fes);
    TestFunction v(fes);

    GridFunction gf(fes);
    auto& x = gf.getData();
    ASSERT_EQ(x.size(), static_cast<Eigen::Index>(fes.getSize()));
    for (Eigen::Index i = 0; i < x.size(); ++i)
      x(i) = static_cast<Real>(i + 1);

    BilinearForm mass(u, v);
    mass = Integral(Dot(u, v));
    mass.assemble();

    LinearForm load(v);
    load = Integral(Dot(gf, v));
    load.assemble();

    const auto residual = mass.getOperator() * x - load.getVector();
    EXPECT_NEAR(residual.norm(), 0.0, 1e-13);
  }
}
