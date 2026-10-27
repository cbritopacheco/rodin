/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Ownership and validation of matrix-valued finite element spaces. */
#include <gtest/gtest.h>
#include <string>
#include <type_traits>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  template <class Scalar>
  class SpaceLifecycle : public ::testing::Test
  {};
  using SpaceScalars = ::testing::Types<Real, Complex>;
  TYPED_TEST_SUITE(SpaceLifecycle, SpaceScalars);

  /// @brief Assignments replace mesh, matrix dimensions and owned maps consistently.
  TYPED_TEST(SpaceLifecycle, MatrixCopiesAndAssignmentsRecoverConstantFields)
  {
    using Scalar = TypeParam;
    using Matrix = Math::SpatialMatrix<Scalar>;
    auto originalMesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {3, 3});
    auto otherMesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    for (auto* mesh : {&originalMesh, &otherMesh})
    {
      mesh->getConnectivity().compute(1, 0);
      mesh->getConnectivity().compute(1, 2);
      mesh->getConnectivity().compute(2, 1);
    }
    Matrix expected(2, 3);
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, -2);
      else
        return Scalar(1);
    }();
    for (size_t row = 0; row < 2; ++row)
    {
      for (size_t column = 0; column < 3; ++column)
        expected(row, column) = factor * Scalar(1 + 3 * row + column);
    }
    auto check = [&](auto factory) {
      using Space = decltype(factory(originalMesh, size_t(2), size_t(3)));
      Space copied = [&] {
        Space original = factory(originalMesh, 2, 3);
        return Space(original);
      }();
      Space moved(std::move(copied));
      Space assigned = factory(otherMesh, 3, 1);
      assigned = moved;
      Space moveAssigned = factory(otherMesh, 1, 1);
      moveAssigned = std::move(assigned);
      auto verify = [&](const Space& space) {
        EXPECT_EQ(&space.getMesh(), &originalMesh);
        EXPECT_EQ(&space.getScalarSpace().getMesh(), &originalMesh);
        EXPECT_EQ(space.getRows(), 2);
        EXPECT_EQ(space.getColumns(), 3);
        EXPECT_EQ(space.getVectorDimension(), 6);
        EXPECT_EQ(space.getSize(), 6 * space.getScalarSpace().getSize());
        GridFunction field(space);
        field = MatrixFunction(expected);
        for (auto cell = originalMesh.getCell(); cell; ++cell)
        {
          const auto& dofs = space.getDOFs(2, cell->getIndex());
          const auto& element = space.getFiniteElement(2, cell->getIndex());
          EXPECT_EQ(static_cast<size_t>(dofs.size()), element.getCount());
          for (Index local = 0; local < static_cast<Index>(dofs.size()); ++local)
            EXPECT_EQ(space.getGlobalIndex({2, cell->getIndex()}, local), dofs[local]);
          const Point point(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
          EXPECT_LT((field(point) - expected).norm(), 1e-11);
        }
      };
      verify(moved);
      verify(moveAssigned);
      moved = factory(otherMesh, 1, 1);
      verify(moveAssigned);
    };
    check([](const LocalMesh& mesh, size_t rows, size_t columns) {
      return P0<Matrix>(mesh, rows, columns);
    });
    check([](const LocalMesh& mesh, size_t rows, size_t columns) {
      return P0g<Matrix, LocalMesh>(mesh, rows, columns);
    });
    check([](const LocalMesh& mesh, size_t rows, size_t columns) {
      return P1<Matrix>(mesh, rows, columns);
    });
    check([](const LocalMesh& mesh, size_t rows, size_t columns) {
      return H1<2, Matrix>(std::integral_constant<size_t, 2>{}, mesh, rows, columns);
    });
  }

  /// @brief Invalid matrix shapes raise native exceptions identifying the owning member.
  TYPED_TEST(SpaceLifecycle, InvalidMatrixShapesHaveMemberContext)
  {
    using Matrix = Math::SpatialMatrix<TypeParam>;
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(1, 0);
    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(2, 1);
    auto check = [&](auto factory) {
      for (const auto& [rows, columns] :
        {std::pair<size_t, size_t>{0, 1}, {4, 1}, {1, 0}, {1, 4}})
      {
        try
        {
          auto invalid = factory(rows, columns);
          FAIL() << "Invalid matrix shape accepted";
        }
        catch (const Alert::Exception& error)
        {
          const std::string message = error.what();
          EXPECT_NE(message.find("In member function "), std::string::npos);
          EXPECT_NE(message.find("of class "), std::string::npos);
          EXPECT_NE(message.find("SpatialMatrix ranges require 1 to 3 rows and columns."),
            std::string::npos);
        }
      }
    };
    check([&](size_t rows, size_t columns) { return P0<Matrix>(mesh, rows, columns); });
    check([&](size_t rows, size_t columns) {
      return P0g<Matrix, LocalMesh>(mesh, rows, columns);
    });
    check([&](size_t rows, size_t columns) { return P1<Matrix>(mesh, rows, columns); });
    check([&](size_t rows, size_t columns) {
      return H1<2, Matrix>(std::integral_constant<size_t, 2>{}, mesh, rows, columns);
    });
  }
}

namespace Rodin::Tests::Unit
{
  /// @brief Runtime-selected spaces expose complete, consistent global and local DOF maps.
  TYPED_TEST(SpaceLifecycle, PolymorphicMetadataAndDofMapsMatchTopology)
  {
    using Scalar = TypeParam;
    using Vector = Math::SpatialVector<Scalar>;
    using Matrix = Math::SpatialMatrix<Scalar>;
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    mesh.getConnectivity().compute(1, 0);
    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(2, 1);
    struct Entry
    {
      std::unique_ptr<FiniteElementSpaceBase> space;
      size_t components;
      size_t scalarSize;
      size_t cellBasisCount;
      bool traces;
    };
    std::vector<Entry> spaces;
    auto add = [&](auto space, size_t components, size_t size, size_t cellBasisCount,
                 bool traces) {
      using Space = std::decay_t<decltype(space)>;
      spaces.push_back({std::make_unique<Space>(std::move(space)), components, size,
        cellBasisCount, traces});
    };
    const size_t cells = mesh.getCellCount();
    const size_t vertices = mesh.getVertexCount();
    const size_t quadraticSize =
      vertices + mesh.getPolytopeCount(Polytope::Type::Segment);
    add(P0<Scalar>(mesh), 1, cells, 1, false);
    add(P0<Vector>(mesh, 2), 2, cells, 1, false);
    add(P0<Matrix>(mesh, 2, 3), 6, cells, 1, false);
    add(P0g<Scalar, LocalMesh>(mesh), 1, 1, 1, true);
    add(P0g<Vector, LocalMesh>(mesh, 2), 2, 1, 1, true);
    add(P0g<Matrix, LocalMesh>(mesh, 2, 3), 6, 1, 1, true);
    add(P1<Scalar>(mesh), 1, vertices, 3, true);
    add(P1<Vector>(mesh, 2), 2, vertices, 3, true);
    add(P1<Matrix>(mesh, 2, 3), 6, vertices, 3, true);
    add(H1<2, Scalar>(std::integral_constant<size_t, 2>{}, mesh), 1, quadraticSize, 6,
      true);
    add(H1<2, Vector>(std::integral_constant<size_t, 2>{}, mesh, 2), 2, quadraticSize, 6,
      true);
    add(H1<2, Matrix>(std::integral_constant<size_t, 2>{}, mesh, 2, 3), 6, quadraticSize,
      6, true);
    for (const auto& entry : spaces)
    {
      const auto& space = *entry.space;
      EXPECT_EQ(&space.getMesh(), static_cast<const MeshBase*>(&mesh));
      EXPECT_EQ(space.getVectorDimension(), entry.components);
      EXPECT_EQ(space.getSize(), entry.scalarSize * entry.components);
      EXPECT_TRUE(space == space);
      FlatSet<Index> global;
      for (size_t dimension = entry.traces ? 0 : 2; dimension <= 2; ++dimension)
      {
        for (auto polytope = mesh.getPolytope(dimension); polytope; ++polytope)
        {
          const auto& dofs = space.getDOFs(dimension, polytope->getIndex());
          if (dimension == 2)
          {
            EXPECT_EQ(
              static_cast<size_t>(dofs.size()), entry.cellBasisCount * entry.components);
          }
          FlatSet<Index> local;
          for (Index index = 0; index < static_cast<Index>(dofs.size()); ++index)
          {
            EXPECT_LT(dofs[index], space.getSize());
            EXPECT_TRUE(local.insert(dofs[index]).second);
            EXPECT_EQ(space.getGlobalIndex({dimension, polytope->getIndex()}, index),
              dofs[index]);
            global.insert(dofs[index]);
          }
        }
      }
      EXPECT_EQ(global.size(), space.getSize());
      for (const auto& other : spaces)
      {
        EXPECT_EQ(space == *other.space, entry.space.get() == other.space.get());
        EXPECT_EQ(space != *other.space, entry.space.get() != other.space.get());
      }
    }
  }
}
