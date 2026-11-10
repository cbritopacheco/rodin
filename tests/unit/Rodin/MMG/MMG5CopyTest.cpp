/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/MMG.h"

using namespace Rodin;
using namespace Rodin::Geometry;

/// @brief Native copies retain the capacity and free lists needed for insertion.
TEST(MMG5Copy, InsertionAfterSourceDestruction)
{
  using Owner = std::unique_ptr<MMG5_Mesh, decltype(&MMG::V5::Mesh::destroy)>;
  for (const bool volume : {false, true})
  {
    SCOPED_TRACE(volume);
    auto mesh = volume ? LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2})
                       : LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    Owner source(MMG::V5::Mesh::toV5(mesh), MMG::V5::Mesh::destroy);
    Owner copy(
      static_cast<MMG5_pMesh>(mycalloc(1, sizeof(MMG5_Mesh))), MMG::V5::Mesh::destroy);
    MMG::V5::Mesh::copy(source.get(), copy.get());
    ASSERT_EQ(copy->npmax, source->npmax);
    ASSERT_GT(copy->npmax, copy->np);
    const auto oldCount = copy->np;
    source.reset();
    double coordinates[3] = {0.25, 0.25, 0.25};
    const auto vertex = volume ? MMG3D_newPt(copy.get(), coordinates, 0, 0)
                               : MMG2D_newPt(copy.get(), coordinates, 0);
    ASSERT_GT(vertex, oldCount);
    EXPECT_DOUBLE_EQ(copy->point[vertex].c[0], 0.25);
    EXPECT_DOUBLE_EQ(copy->point[vertex].c[1], 0.25);
    if (volume)
    {
      EXPECT_GT(MMG3D_newElt(copy.get()), mesh.getCellCount());
    }
    else
    {
      EXPECT_GT(MMG2D_newElt(copy.get()), mesh.getCellCount());
    }
  }
}

/// @brief Native scalar/vector transfers preserve interleaved nodal coefficients in both directions.
TEST(MMG5Copy, SolutionGridFunctionRoundtrip)
{
  const auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});

  std::unique_ptr<MMG5_Mesh, decltype(&MMG::V5::Mesh::destroy)> nativeMesh(
    MMG::V5::Mesh::toV5(mesh), MMG::V5::Mesh::destroy);
  const auto check = [&](auto& source, auto& destination, size_t components) {
    for (Index i = 0; i < source.getData().size(); ++i)
      source.getData()(i) = Real(i + 1) / 8;
    std::unique_ptr<MMG5_Sol, decltype(&MMG::V5::Solution::destroy)> native(
      MMG::V5::Solution::create(nativeMesh.get(), components),
      MMG::V5::Solution::destroy);
    MMG::V5::Solution::copy(source, native.get());
    ASSERT_EQ(native->np, mesh.getVertexCount());
    for (size_t i = 0; i < mesh.getVertexCount() * components; ++i)
      EXPECT_EQ(native->m[components + i], source.getData()(i));
    MMG::V5::Solution::copy(native.get(), destination);
    EXPECT_EQ(destination.getData(), source.getData());
  };
  Variational::P1 scalarSpace(mesh);
  Variational::P1 vectorSpace(mesh, 2);
  MMG::RealGridFunction scalarSource(scalarSpace), scalarDestination(scalarSpace);
  MMG::VectorGridFunction vectorSource(vectorSpace), vectorDestination(vectorSpace);
  check(scalarSource, scalarDestination, 1);
  check(vectorSource, vectorDestination, 2);
}
