#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

TEST(Rodin_Variational_Scalar_P1, SanityTest_2D_Square_Build)
{
  constexpr size_t vdim = 1;
  constexpr size_t sdim = 2;
  constexpr size_t mdim = 2;

  Mesh mesh =
    Mesh<Rodin::Context::Serial>::Builder()
    .initialize(sdim)
    .nodes(4)
    .vertex({0, 0})
    .vertex({1, 0})
    .vertex({0, 1})
    .vertex({1, 1})
    .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
    .polytope(Polytope::Geometry::Triangle, {1, 3, 2})
    .finalize();

  EXPECT_EQ(mesh.getDimension(), mdim);
  EXPECT_EQ(mesh.getSpaceDimension(), sdim);

  P1 fes(mesh);

  EXPECT_EQ(fes.getVectorDimension(), vdim);
  EXPECT_EQ(fes.getSize(), mesh.getVertexCount());

  EXPECT_EQ(fes.getSize(), 4);
  EXPECT_EQ(fes.getFiniteElement(mdim, 0).getGeometry(), Polytope::Geometry::Triangle);
  EXPECT_EQ(fes.getFiniteElement(mdim, 1).getGeometry(), Polytope::Geometry::Triangle);
}

TEST(Rodin_Variational_Scalar_P1_GridFunction, SanityTest_2D_Square_Project_Sum)
{
  constexpr size_t mdim = 2;

  Mesh mesh =
    Mesh<Rodin::Context::Serial>::Builder()
    .initialize(mdim)
    .nodes(4)
    .vertex({0, 0})
    .vertex({1, 0})
    .vertex({0, 1})
    .vertex({1, 1})
    .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
    .polytope(Polytope::Geometry::Triangle, {1, 3, 2})
    .finalize();

  P1 fes(mesh);
  GridFunction gf(fes);

  ScalarFunction c([](const Geometry::Point& p) { return p.x() + p.y(); } );
  gf.project(c);

  EXPECT_EQ(gf.getRangeShape(), RangeShape(1, 1));

  EXPECT_NEAR(gf.getValue(0), 0, RODIN_FUZZY_CONSTANT);
  EXPECT_NEAR(gf.getValue(1), 1, RODIN_FUZZY_CONSTANT);
  EXPECT_NEAR(gf.getValue(2), 1, RODIN_FUZZY_CONSTANT);
  EXPECT_NEAR(gf.getValue(3), 2, RODIN_FUZZY_CONSTANT);
}

TEST(Rodin_Variational_Scalar_P1_GridFunction_FuzzyTest, TriangularUniformGrid16_ProjectOnBoundary_Constant)
{
  constexpr size_t mdim = 2;

  Mesh mesh = SerialMesh::UniformGrid(Polytope::Geometry::Triangle, 16, 16);

  mesh.getConnectivity().compute(mdim - 1, mdim);

  P1 fes(mesh);
  GridFunction gf(fes);

  ScalarFunction c = 1.0;
  gf.projectOnBoundary(c);

  EXPECT_EQ(gf.getRangeShape(), RangeShape(1, 1));
}

TEST(Rodin_Variational_Scalar_P1_GridFunction, FuzzyTest_TriangularUniformGrid16_ProjectOnBoundary_Sum)
{
  constexpr size_t mdim = 2;

  Mesh mesh = SerialMesh::UniformGrid(Polytope::Geometry::Triangle, 16, 16);

  mesh.getConnectivity().compute(mdim - 1, mdim);

  P1 fes(mesh);
  GridFunction gf(fes);

  ScalarFunction c([](const Geometry::Point& p) { return p.x() + p.y(); } );
  gf.projectOnBoundary(c);

  EXPECT_EQ(gf.getRangeShape(), RangeShape(1, 1));
}

TEST(Rodin_Variational_Scalar_P1_GridFunction, FuzzyTest_2D_Square_Project_LinearFunction)
{
  constexpr size_t mdim = 2;

  Mesh mesh =
    Mesh<Rodin::Context::Serial>::Builder()
    .initialize(mdim)
    .nodes(4)
    .vertex({0, 0})
    .vertex({1, 0})
    .vertex({0, 1})
    .vertex({1, 1})
    .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
    .polytope(Polytope::Geometry::Triangle, {1, 3, 2})
    .finalize();

  P1 fes(mesh);
  GridFunction gf1(fes);

  ScalarFunction c([](const Geometry::Point& p) { return p.x() + p.y(); } );
  gf1.project(c);

  EXPECT_EQ(gf1.getRangeShape(), RangeShape(1, 1));

  RandomFloat gen(0.0, 1.0);
  {
    Index idx = 0;
    auto it = mesh.getPolytope(mdim, idx);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mdim, idx);

    for (size_t i = 0; i < 25; i++)
    {
      const Scalar x = gen();
      const Scalar y = gen();
      const Math::Vector rc{{x, y}};
      const Math::Vector pc = trans.transform(rc);
      const Point p(polytope, trans, Math::Vector{{x, y}});
      EXPECT_NEAR(gf1.getValue(p), pc.x() + pc.y(), RODIN_FUZZY_CONSTANT);
    }
  }

  {
    Index idx = 1;
    auto it = mesh.getPolytope(mdim, idx);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mdim, idx);

    for (size_t i = 0; i < 25; i++)
    {
      const Scalar x = gen();
      const Scalar y = gen();
      const Math::Vector rc{{x, y}};
      const Math::Vector pc = trans.transform(rc);
      const Point p(polytope, trans, rc);
      EXPECT_NEAR(gf1.getValue(p), pc.x() + pc.y(), RODIN_FUZZY_CONSTANT);
    }
  }

  GridFunction gf2(fes);
  gf2.project([](const Geometry::Point& p) { return 5 * p.x() + 100 * p.y(); });

  EXPECT_EQ(gf1.getRangeShape(), RangeShape(1, 1));

  {
    auto it = mesh.getPolytope(mdim, 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mdim, 0);

    for (size_t i = 0; i < 25; i++)
    {
      const Scalar x = gen();
      const Scalar y = gen();
      const Math::Vector rc{{x, y}};
      const Math::Vector pc = trans.transform(rc);
      const Point p(polytope, trans, rc);
      EXPECT_NEAR(gf2.getValue(p), 5 * pc.x() + 100 * pc.y(), RODIN_FUZZY_CONSTANT);
    }
  }

  GridFunction gf3(fes);
  gf3.project([](const Geometry::Point& p) { return p.x() - p.y(); });

  EXPECT_EQ(gf1.getRangeShape(), RangeShape(1, 1));

  {
    auto it = mesh.getPolytope(mdim, 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mdim, 0);

    for (size_t i = 0; i < 25; i++)
    {
      const Scalar x = gen();
      const Scalar y = gen();
      const Math::Vector rc{{x, y}};
      const Math::Vector pc = trans.transform(rc);
      const Point p(polytope, trans, rc);
      EXPECT_NEAR(gf3.getValue(p), pc.x() - pc.y(), RODIN_FUZZY_CONSTANT);
    }
  }

  GridFunction gf4(fes);
  gf4.project([](const Geometry::Point& p) { return 666 * p.x() - 999 * p.y(); });

  EXPECT_EQ(gf1.getRangeShape(), RangeShape(1, 1));

  {
    auto it = mesh.getPolytope(mdim, 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mdim, 0);

    for (size_t i = 0; i < 25; i++)
    {
      const Scalar x = gen();
      const Scalar y = gen();
      const Math::Vector rc{{x, y}};
      const Math::Vector pc = trans.transform(rc);
      const Point p(polytope, trans, rc);
      EXPECT_NEAR(gf4.getValue(p), 666 * pc.x() - 999 * pc.y(), RODIN_FUZZY_CONSTANT);
    }
  }
}

TEST(Rodin_Variational_Vector_P1, SanityTest_2D_Square_Build)
{
  constexpr size_t sdim = 2;
  constexpr size_t mdim = sdim;
  constexpr size_t vdim = mdim;

  Mesh mesh =
    Mesh<Rodin::Context::Serial>::Builder()
    .initialize(sdim)
    .nodes(4)
    .vertex({0, 0})
    .vertex({1, 0})
    .vertex({0, 1})
    .vertex({1, 1})
    .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
    .polytope(Polytope::Geometry::Triangle, {1, 3, 2})
    .finalize();

  EXPECT_EQ(mesh.getDimension(), mdim);
  EXPECT_EQ(mesh.getSpaceDimension(), sdim);

  P1 fes(mesh, mdim);

  EXPECT_EQ(fes.getVectorDimension(), vdim);
  EXPECT_EQ(fes.getSize(), vdim * mesh.getVertexCount());
  EXPECT_EQ(fes.getFiniteElement(mdim, 0).getGeometry(), Polytope::Geometry::Triangle);
  EXPECT_EQ(fes.getFiniteElement(mdim, 1).getGeometry(), Polytope::Geometry::Triangle);
}

TEST(Rodin_Variational_Vector_P1_GridFunction, FuzzyTest_2D_Square_Project)
{
  constexpr size_t mdim = 2;

  Mesh mesh =
    Mesh<Rodin::Context::Serial>::Builder()
    .initialize(mdim)
    .nodes(4)
    .vertex({0, 0})
    .vertex({1, 0})
    .vertex({0, 1})
    .vertex({1, 1})
    .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
    .polytope(Polytope::Geometry::Triangle, {1, 3, 2})
    .finalize();

  P1 fes(mesh, mdim);
  GridFunction gf1(fes);

  EXPECT_EQ(gf1.getRangeShape(), RangeShape(mdim, 1));

  VectorFunction c1 = {
    [](const Geometry::Point& p){ return p.x(); },
    [](const Geometry::Point& p){ return p.y(); }
  };

  gf1.project(c1);

  RandomFloat gen(0.0, 1.0);

  {
    auto it = mesh.getPolytope(mdim, 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mdim, 0);

    for (size_t i = 0; i < 25; i++)
    {
      const Scalar x = gen();
      const Scalar y = gen();
      const Math::Vector rc{{x, y}};
      const Math::Vector pc = trans.transform(rc);
      const Point p(polytope, trans, rc);
      EXPECT_NEAR((gf1.getValue(p) - Math::Vector{{pc.x(), pc.y()}}).norm(), 0, RODIN_FUZZY_CONSTANT);
      break;
    }
  }

  GridFunction gf2(fes);

  EXPECT_EQ(gf2.getRangeShape(), RangeShape(mdim, 1));

  VectorFunction c2 = {
    [](const Geometry::Point& p){ return p.y(); },
    [](const Geometry::Point& p){ return p.x(); }
  };

  gf2.project(c2);

  {
    auto it = mesh.getPolytope(mdim, 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mdim, 0);

    for (size_t i = 0; i < 25; i++)
    {
      const Scalar x = gen();
      const Scalar y = gen();
      const Math::Vector rc{{x, y}};
      const Math::Vector pc = trans.transform(rc);
      const Point p(polytope, trans, rc);
      EXPECT_NEAR((gf2.getValue(p) - Math::Vector{{pc.y(), pc.x()}}).norm(), 0, RODIN_FUZZY_CONSTANT);
    }
  }

  GridFunction gf3(fes);

  EXPECT_EQ(gf3.getRangeShape(), RangeShape(mdim, 1));

  VectorFunction c3 = { [](const Geometry::Point& p){ return p.x() + p.y(); }, 0 };

  gf3.project(c3);

  {
    auto it = mesh.getPolytope(mdim, 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mdim, 0);

    for (size_t i = 0; i < 25; i++)
    {
      const Scalar x = gen();
      const Scalar y = gen();
      const Math::Vector rc{{x, y}};
      const Math::Vector pc = trans.transform(rc);
      const Point p(polytope, trans, rc);
      EXPECT_NEAR((gf3.getValue(p) - Math::Vector{{pc.x() + pc.y(), 0}}).norm(), 0, RODIN_FUZZY_CONSTANT);
    }
  }

  GridFunction gf4(fes);

  EXPECT_EQ(gf4.getRangeShape(), RangeShape(mdim, 1));

  VectorFunction c4 = {
    [](const Geometry::Point& p){ return 999 * p.x() - 100 * p.y(); },
    [](const Geometry::Point& p){ return -5 * p.x() + 666 * p.y(); },
  };

  gf4.project(c4);

  {
    auto it = mesh.getPolytope(mdim, 0);
    const auto& polytope = *it;
    const auto& trans = mesh.getPolytopeTransformation(mdim, 0);

    for (size_t i = 0; i < 25; i++)
    {
      const Scalar x = gen();
      const Scalar y = gen();
      const Math::Vector rc{{x, y}};
      const Math::Vector pc = trans.transform(rc);
      const Point p(polytope, trans, rc);
      const Math::Vector actual{{999 * p.x() - 100 * p.y(), -5 * p.x() + 666 * p.y()}};
      EXPECT_NEAR((gf4.getValue(p) - actual).norm(), 0, RODIN_FUZZY_CONSTANT);
    }
  }
}

TEST(Rodin_Variational_Scalar_P1_TrialFunction, FuzzyTest_UniformGrid_4x4)
{
  Mesh mesh = SerialMesh::UniformGrid(Polytope::Geometry::Triangle, 4, 4);
  P1 fes(mesh);
  TrialFunction u(fes);

  EXPECT_EQ(u.getRangeShape(), RangeShape(1, 1));
}

TEST(Rodin_Variational_Scalar_P1_TestFunction, FuzzyTest_UniformGrid_4x4)
{
  Mesh mesh = SerialMesh::UniformGrid(Polytope::Geometry::Triangle, 4, 4);
  P1 fes(mesh);
  TrialFunction v(fes);

  EXPECT_EQ(v.getRangeShape(), RangeShape(1, 1));
}


TEST(Rodin_Variational_Scalar_P1_LinearForm, FuzzyTest_UniformGrid_4x4)
{
  Mesh mesh = SerialMesh::UniformGrid(Polytope::Geometry::Triangle, 4, 4);
  P1 fes(mesh);
  TestFunction v(fes);
  LinearForm lf(v);
  lf = Integral(v);
  lf.assemble();
}

