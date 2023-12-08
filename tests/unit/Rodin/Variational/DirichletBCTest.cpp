#include <gtest/gtest.h>

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Rodin_Variational_Scalar_P1_SanityTest, TriangularUniformGrid2)
{
  Mesh mesh = SequentialMesh::UniformGrid(Polytope::Type::Triangle, 2, 2);
  const size_t D = mesh.getDimension();

  mesh.getConnectivity().compute(D - 1, D);

  P1 fes(mesh);

  TrialFunction u(fes);
  u.emplace();

  ScalarFunction c = 1;

  DirichletBC dbc(u, c);
  dbc.assemble();

  EXPECT_EQ(dbc.getDOFs().size(), 4);
}

TEST(Rodin_Variational_Scalar_P1_SanityTest, TriangularUniformGrid16)
{
  Mesh mesh = SequentialMesh::UniformGrid(Polytope::Type::Triangle, 16, 16);
  const size_t D = mesh.getDimension();
  const Attribute attr = RODIN_DEFAULT_POLYTOPE_ATTRIBUTE;

  mesh.getConnectivity().compute(D - 1, D);

  P1 fes(mesh);

  TrialFunction u(fes);
  u.emplace();

  ScalarFunction c = 1;

  DirichletBC dbc(u, c);
  dbc.on(attr);
  dbc.assemble();

  EXPECT_EQ(dbc.getDOFs().size(), 256);
}
