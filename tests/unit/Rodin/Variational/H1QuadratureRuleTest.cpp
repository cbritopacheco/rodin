#include <gtest/gtest.h>

#include "Rodin/Variational.h"
#include "Rodin/Assembly/Default.h"
#include "Rodin/Geometry/Mesh.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
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

  TEST(Rodin_Variational_H1QuadratureRule, MixedOrder_Jacobian_Assembles)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 0);
    const size_t vdim = 2;
    H1<2, Math::Vector<Real>> fesTr(std::integral_constant<size_t, 2>{}, mesh, vdim);
    H1<1, Math::Vector<Real>> fesTe(std::integral_constant<size_t, 1>{}, mesh, vdim);

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
}
