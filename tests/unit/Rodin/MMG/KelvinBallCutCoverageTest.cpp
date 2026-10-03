/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "examples/ShapeOptimization/KelvinBall/RotatedNitscheIntegrator.h"

namespace KelvinBall
{
  class CutCoverageTest : public ::testing::Test
  {
    protected:
      Mesh makeMesh(Real shift)
      {
        Mesh::Builder builder;
        builder.initialize(3).nodes(8);
        for (const auto& x : std::array<std::array<Real, 3>, 8>{{
          {0, 0, 0}, {1, 0, 0}, {0, 1, -1}, {0.2, 0.5, 0},
          {shift, 0, 0}, {shift + 1, 0, 0}, {shift, 1, 1},
          {shift + 0.2, 0.5, 0}}})
        {
          Math::SpatialPoint point(3);
          for (size_t i = 0; i < 3; ++i)
            point(i) = x[i];
          builder.vertex(point);
        }
        Index cell;
        IndexArray vertices(4);
        vertices << 0, 1, 2, 3;
        builder.polytope(Polytope::Type::Tetrahedron, vertices, cell);
        builder.attribute({3, cell}, Fluid);
        vertices << 5, 4, 6, 7;
        builder.polytope(Polytope::Type::Tetrahedron, vertices, cell);
        builder.attribute({3, cell}, Fluid);
        Mesh mesh = builder.finalize();
        mesh.getConnectivity().compute(2, 3);
        for (auto face = mesh.getBoundary(); face; ++face)
        {
          bool slave = true, master = true;
          for (Index v : face->getVertices())
          {
            const auto& x = mesh.getVertexCoordinates(v);
            slave = slave && std::abs(x(1) + x(2)) < 1e-12;
            master = master && std::abs(x(1) - x(2)) < 1e-12;
          }
          mesh.setAttribute({2, face->getIndex()},
            slave ? SigmaMinus : master ? SigmaPlus : Outer);
        }
        prepare(mesh);
        return mesh;
      }

      using System = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

      System makeSystem(size_t size)
      {
        System system;
        system.getOperator().resize(size, size);
        system.getVector().setZero(size);
        system.getSolution().setZero(size);
        return system;
      }
  };

  TEST_F(CutCoverageTest, CompletelyMissingTraceAddsNeitherMatrixNorLoad)
  {
    auto mesh = makeMesh(2);
    P1 scalar(mesh);
    P1 vector(mesh, 3);
    GridFunction reference(scalar);
    reference = RealFunction([](const Geometry::Point& p) {
      return p.getPhysicalCoordinates()(0);
    });
    RotatedNitscheIntegrator coupling(mesh);
    auto scalarSystem = makeSystem(scalar.getSize());
    EXPECT_NO_THROW(coupling.assembleScalarTracePenalty(
      scalar, scalarSystem, 1, reference));
    EXPECT_EQ(scalarSystem.getOperator().norm(), 0);
    EXPECT_EQ(scalarSystem.getVector().norm(), 0);
    auto vectorSystem = makeSystem(vector.getSize());
    EXPECT_NO_THROW(coupling.assembleVector(vector, vectorSystem, 1, 10, {}));
    EXPECT_EQ(vectorSystem.getOperator().norm(), 0);
    const size_t block = vector.getSize() + scalar.getSize();
    const std::array<size_t, 6> offsets{0, vector.getSize(), block,
      block + vector.getSize(), 2 * block, 2 * block + vector.getSize()};
    auto stokesSystem = makeSystem(3 * block);
    EXPECT_NO_THROW(coupling.assembleStokes<1>(
      vector, scalar, offsets, stokesSystem, 1, 10, 0.05, {}));
    EXPECT_EQ(stokesSystem.getOperator().norm(), 0);
    EXPECT_EQ(stokesSystem.getVector().norm(), 0);
    GridFunction v(vector);
    v = VectorFunction(0.0, 0.0, 0.0);
    EXPECT_NO_THROW(coupling.scalarJump(reference));
    EXPECT_NO_THROW(coupling.vectorJump(v));
    EXPECT_NO_THROW(coupling.familyJump(v, v, v));
  }

  TEST_F(CutCoverageTest, PartialOverlapRetainsSymmetricNonzeroCoupling)
  {
    auto mesh = makeMesh(0.3);
    P1 scalar(mesh);
    RotatedNitscheIntegrator coupling(mesh);
    size_t matched = 0, missed = 0;
    for (const auto& pair : RotationPairs)
    {
      for (auto face = mesh.getBoundary(); face; ++face)
      {
        if (face->getAttribute() != pair.slave)
          continue;
        const auto& formula = QF::PolytopeQuadratureFormula::get(4, face->getGeometry());
        const auto& quadrature = face->getQuadrature(formula);
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          if (coupling.getLocator().locate(
                pair.master, pair.rotation * quadrature.getPoint(qp).vector()))
            ++matched;
          else
            ++missed;
        }
      }
    }
    ASSERT_GT(matched, 0);
    ASSERT_GT(missed, 0);
    auto system = makeSystem(scalar.getSize());
    EXPECT_NO_THROW(coupling.assembleScalarTracePenalty(
      scalar, system, 1, FlatSet<Attribute>{}));
    EXPECT_GT(system.getOperator().norm(), 0);
    EXPECT_LT((system.getOperator() -
      Math::SparseMatrix<Real>(system.getOperator().transpose())).norm(), 1e-12);
    Math::Vector<Real> ones = Math::Vector<Real>::Ones(scalar.getSize());
    EXPECT_LT((system.getOperator() * ones).norm(), 1e-12);
  }
}
