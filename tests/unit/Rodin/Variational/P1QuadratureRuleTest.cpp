#include <gtest/gtest.h>

#include "Rodin/Variational.h"
#include "Rodin/Assembly.h"
#include "Rodin/Geometry/Mesh.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_P1QuadratureRule, MixedSpaces_GradGrad_Assembles)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 fesTr(mesh);
    P1 fesTe(mesh); // distinct space instance to trigger mixed-space path

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

  TEST(Rodin_Variational_P1QuadratureRule, MixedSpaces_VectorMass_Assembles)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const size_t vdim = 2;
    P1<Math::Vector<Real>> fesTr(mesh, vdim);
    P1<Math::Vector<Real>> fesTe(mesh, vdim);

    TrialFunction u(fesTr);
    TestFunction v(fesTe);

    RealFunction coeff(2.0);

    BilinearForm bf(u, v);
    bf = Integral(Dot(Mult(coeff, u), v));
    bf.assemble();

    const auto& mat = bf.getOperator();
    EXPECT_EQ(mat.rows(), static_cast<Eigen::Index>(fesTe.getSize()));
    EXPECT_EQ(mat.cols(), static_cast<Eigen::Index>(fesTr.getSize()));
    EXPECT_GT(mat.norm(), 0.0);
  }
}
