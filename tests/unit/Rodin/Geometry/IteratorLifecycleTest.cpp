/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Iterator ownership and typed-assignment contracts. */
#include <gtest/gtest.h>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>
#include "Rodin/Geometry.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  Mesh<Context::Local> referenceMesh(Polytope::Type geometry)
  {
    const Polytope::Traits traits(geometry);
    Mesh<Context::Local>::Builder builder;
    builder.initialize(traits.getDimension()).nodes(traits.getVertexCount());
    IndexArray vertices(traits.getVertexCount());
    for (size_t a = 0; a < traits.getVertexCount(); ++a)
    {
      builder.vertex(traits.getVertex(a));
      vertices(a) = a;
    }
    builder.polytope(geometry, vertices);
    auto mesh = builder.finalize();
    for (size_t k = traits.getDimension(); k > 0; --k)
      mesh.getConnectivity().compute(k, k - 1);
    for (size_t k = 0; k <= traits.getDimension(); ++k)
    {
      for (Index i = 0; i < mesh.getPolytopeCount(k); ++i)
        mesh.setAttribute({k, i}, 100 * k + i);
    }
    return mesh;
  }

  void checkEntity(
    const Polytope& entity, const MeshBase& mesh, size_t dimension, Index index)
  {
    const auto expected = mesh.getPolytope(dimension, index);
    EXPECT_EQ(&entity.getMesh(), &mesh);
    EXPECT_EQ(entity.getDimension(), dimension);
    EXPECT_EQ(entity.getIndex(), index);
    EXPECT_EQ(entity.getGeometry(), expected->getGeometry());
    EXPECT_EQ(entity.getAttribute(), Optional<Attribute>(100 * dimension + index));
    const auto& actualVertices = entity.getVertices();
    const auto& expectedVertices = expected->getVertices();
    EXPECT_EQ(std::vector<Index>(actualVertices.begin(), actualVertices.end()),
      std::vector<Index>(expectedVertices.begin(), expectedVertices.end()));
  }

  template <class Factory>
  void checkReleasedOwnership(const MeshBase& mesh, size_t dimension, Factory make)
  {
    using Iterator = decltype(make(std::vector<Index>{}));
    using Entity = std::remove_pointer_t<decltype(std::declval<Iterator&>().release())>;
    const Index last = mesh.getPolytopeCount(dimension) - 1;
    const std::vector<Index> order{last, 0, last};
    std::vector<std::unique_ptr<Entity>> released;
    {
      auto iterator = make(std::vector<Index>(order));
      // Cold release generates an independently owned entity without advancing.
      released.emplace_back(iterator.release());
      ASSERT_TRUE(released.back());
      EXPECT_EQ(iterator.key().second, last);
      checkEntity(*iterator, mesh, dimension, last);
      // Warm release transfers the cached value, and dereferencing recreates it.
      released.emplace_back(iterator.release());
      ASSERT_TRUE(released.back());
      checkEntity(*iterator, mesh, dimension, last);
      ++iterator;
      // Dirty release must generate the new entity rather than the cached one.
      released.emplace_back(iterator.release());
      ASSERT_TRUE(released.back());
      checkEntity(*iterator, mesh, dimension, 0);
      ++iterator;
      checkEntity(*iterator, mesh, dimension, last);
      ++iterator;
      EXPECT_TRUE(iterator.end());
      EXPECT_EQ(iterator.release(), nullptr);
    }
    // All released entities remain usable after their iterator is destroyed.
    ASSERT_EQ(released.size(), 3);
    checkEntity(*released[0], mesh, dimension, last);
    checkEntity(*released[1], mesh, dimension, last);
    checkEntity(*released[2], mesh, dimension, 0);
    auto empty = make(std::vector<Index>{});
    EXPECT_TRUE(empty.end());
    EXPECT_EQ(empty.release(), nullptr);
  }

  template <class Iterator>
  void checkTypedAssignment(const MeshBase& mesh, size_t dimension)
  {
    const Index last = mesh.getPolytopeCount(dimension) - 1;
    PolytopeIterator target(0, mesh, VectorIndexGenerator({0}));
    checkEntity(*target, mesh, 0, 0);
    {
      Iterator source(mesh, VectorIndexGenerator({0, last, 0}));
      checkEntity(*source, mesh, dimension, 0);
      ++source;
      target = std::move(source);
    }
    EXPECT_EQ(target.key().first, dimension);
    EXPECT_EQ(target.key().second, last);
    checkEntity(*target, mesh, dimension, last);
    ++target;
    checkEntity(*target, mesh, dimension, 0);
    ++target;
    EXPECT_TRUE(target.end());
    EXPECT_EQ(target.release(), nullptr);
  }
}

/// @brief Cold, cached and dirty releases own values independently for every entity kind.
TEST(Geometry_IteratorLifecycle, ReleasedEntitiesSurviveIteratorDestruction)
{
  for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(int(geometry));
    const auto mesh = referenceMesh(geometry);
    for (size_t k = 0; k <= mesh.getDimension(); ++k)
    {
      checkReleasedOwnership(mesh, k, [&](std::vector<Index> indices) {
        return PolytopeIterator(k, mesh, VectorIndexGenerator(std::move(indices)));
      });
    }
    checkReleasedOwnership(mesh, mesh.getDimension(), [&](std::vector<Index> indices) {
      return CellIterator(mesh, VectorIndexGenerator(std::move(indices)));
    });
    checkReleasedOwnership(
      mesh, mesh.getDimension() - 1, [&](std::vector<Index> indices) {
        return FaceIterator(mesh, VectorIndexGenerator(std::move(indices)));
      });
    checkReleasedOwnership(mesh, 0, [&](std::vector<Index> indices) {
      return VertexIterator(mesh, VectorIndexGenerator(std::move(indices)));
    });
  }
}

/// @brief Typed assignment preserves advanced cursors and replaces a warm generic cache.
TEST(Geometry_IteratorLifecycle, TypedAssignmentPreservesCursorAndInvalidatesCache)
{
  for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(int(geometry));
    const auto mesh = referenceMesh(geometry);
    checkTypedAssignment<CellIterator>(mesh, mesh.getDimension());
    checkTypedAssignment<FaceIterator>(mesh, mesh.getDimension() - 1);
    checkTypedAssignment<VertexIterator>(mesh, 0);
  }
}
