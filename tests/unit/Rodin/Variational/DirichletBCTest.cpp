#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_Real_P1_SanityTest, TriangularUniformGrid2)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    const size_t D = mesh.getDimension();

    mesh.getConnectivity().compute(D - 1, D);

    P1 fes(mesh);

    TrialFunction u(fes);

    RealFunction c = 1;

    DirichletBC dbc(u, c);
    dbc.assemble();

    EXPECT_EQ(dbc.getDOFs().size(), 4);
  }

  TEST(Rodin_Variational_Real_P1_SanityTest, TriangularUniformGrid16)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 16, 16 });
    const size_t D = mesh.getDimension();
    const Attribute attr = 1;

    mesh.getConnectivity().compute(D - 1, D);

    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      mesh.setAttribute(it.key(), attr);
    }

    P1 fes(mesh);

    TrialFunction u(fes);

    RealFunction c = 1;

    DirichletBC dbc(u, c);
    dbc.on(attr);
    dbc.assemble();

    EXPECT_EQ(dbc.getDOFs().size(), 60);
  }

  TEST(Rodin_Variational_ShapeFunctionSetPoint, H1MatchesBasis)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D - 1, 0);

    H1<2, Real> fes(std::integral_constant<size_t, 2>{}, mesh);
    TestFunction v(fes);

    Math::SpatialVector<Real> rc(2);
    rc[0] = 0.2;
    rc[1] = 0.3;

    auto it = mesh.getPolytope(D, 0);
    const Geometry::Point p(*it, rc);

    v.setPoint(p);

    const auto& fe = fes.getFiniteElement(D, 0);
    for (size_t local = 0; local < fe.getCount(); local++)
      EXPECT_NEAR(v.getBasis(local), fe.getBasis(local)(rc), 1e-14);
  }

  TEST(Rodin_Variational_ShapeFunctionSetPoint, H1GradMatchesBasisGradient)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D - 1, 0);

    H1<2, Real> fes(std::integral_constant<size_t, 2>{}, mesh);
    TestFunction v(fes);
    auto gv = Grad(v);

    Math::SpatialVector<Real> rc(2);
    rc[0] = 0.2;
    rc[1] = 0.3;

    auto it = mesh.getPolytope(D, 0);
    const Geometry::Point p(*it, rc);

    gv.setPoint(p);

    const auto& fe = fes.getFiniteElement(D, 0);
    const auto JinvT = p.getJacobianInverse().transpose();
    for (size_t local = 0; local < fe.getCount(); local++)
    {
      Math::SpatialVector<Real> ref(D);
      for (size_t d = 0; d < D; d++)
        ref(d) = fe.getBasis(local).template getDerivative<1>(d)(rc);

      const auto expected = JinvT * ref;
      const auto actual = gv.getBasis(local);
      for (size_t d = 0; d < D; d++)
        EXPECT_NEAR(actual(d), expected(d), 1e-14);
    }
  }
}
