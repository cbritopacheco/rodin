/*
 * Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or copy at https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <boost/mpi/collectives.hpp>
#include <Rodin/MPI/Location.h>
#include <Rodin/MPI/Geometry/SubMesh.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/Variational.h>
#include <Rodin/MPI/Variational/P0.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  boost::mpi::environment* environment;
  boost::mpi::communicator* world;

  MPIMesh overlappingMesh(const Context::MPI& context)
  {
    Shard::Builder builder;
    builder.initialize(1, 2);
    // Distinct topological entities with overlapping physical images.
    for (Index c = 0; c < 3; ++c)
    {
      builder.vertex(2 * c, Math::SpatialPoint::Zero(2), Shard::State::Shared);
      Math::SpatialPoint end(2);
      end[0] = 1;
      end[1] = 0;
      builder.vertex(2 * c + 1, end, Shard::State::Ghost);
      builder.setOwner(0, 2 * c, 0).setOwner(0, 2 * c + 1, 0);
      IndexArray vertices(2);
      vertices << 2 * c, 2 * c + 1;
      const auto state = c == 0 ? Shard::State::Ghost
        : c == 1                ? Shard::State::Shared
                                : Shard::State::Owned;
      builder.polytope(1, 100 + c, Polytope::Type::Segment, vertices, state);
      if (state != Shard::State::Owned)
        builder.setOwner(1, c, 0);
    }
    return MPIMesh::Builder(context).initialize(builder.finalize()).finalize();
  }

  // Each cell is disjoint and uses the reference geometry. Extra ambient axes
  // exercise embedded transformations. Compression and curvature stay moderate
  // so these fixtures isolate ownership/transfer rather than inverse conditioning.
  template <size_t K>
  LocalMesh curvedMesh(Polytope::Type type, bool embedded)
  {
    const Polytope::Traits traits(type);
    const size_t d = traits.getDimension();
    const size_t sdim = embedded ? std::min(size_t(3), d + 1) : std::max(size_t(1), d);
    constexpr Real CellSpacing = Real(2);
    constexpr Real AxisCompression = Real(0.4);
    constexpr Real Curvature = Real(0.12);
    LocalMesh::Builder builder;
    builder.initialize(sdim);
    const size_t cells = static_cast<size_t>(world->size());
    const size_t nv = traits.getVertexCount();
    builder.nodes(cells * nv);
    for (size_t c = 0; c < cells; ++c)
    {
      IndexArray vertices(nv);
      for (size_t v = 0; v < nv; ++v)
      {
        Math::SpatialPoint x = Math::SpatialPoint::Zero(sdim);
        if (d > 0)
          for (size_t j = 0; j < d; ++j)
            x[j] = traits.getVertex(v)[j];
        x[0] = CellSpacing * c + AxisCompression * x[0];
        builder.vertex(x);
        vertices[v] = c * nv + v;
      }
      if (d > 0)
        builder.polytope(type, vertices);
    }
    auto mesh = builder.finalize();
    if (d > 0)
    {
      Variational::RealH1Element<K> element(type);
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        PointCloud nodes(sdim, element.getCount());
        for (size_t a = 0; a < element.getCount(); ++a)
        {
          Math::SpatialPoint x;
          cell->getTransformation().transform(x, element.getNode(a));
          const Real t = element.getNode(a)[0];
          if constexpr (K > 1)
            x[sdim - 1] += Curvature * t * (1 - t);
          for (size_t j = 0; j < sdim; ++j)
            nodes(j, a) = x[j];
        }
        mesh.setPolytopeTransformation({d, cell->getIndex()},
          new ParametricTransformation<Variational::RealH1Element<K>>(
            std::move(nodes), element));
      }
      mesh.getConnectivity().compute(d, d);
    }
    return mesh;
  }

  template <size_t K>
  void checkDistributedGeometry(
    Polytope::Type type, bool embedded, bool defaultMaps = false)
  {
    Context::MPI context(*environment, *world);
    Sharder<Context::MPI> sharder(context);
    const size_t d = Polytope::Traits(type).getDimension();
    Math::SpatialPoint expected;
    std::vector<Real> expectedJacobian;
    Math::SpatialPoint rc = d == 0 ? Math::SpatialPoint(0)
                                   : (Real(0.75) * Polytope::Traits(type).getCentroid() +
                                       Real(0.25) * Polytope::Traits(type).getVertex(0));
    if (world->rank() == 0)
    {
      auto parent = curvedMesh<K>(type, embedded);
      if (defaultMaps)
      {
        parent.flush();
        // Cache every P1 chart, including points, so remote shard transport
        // must serialize these polymorphic maps rather than rebuild them.
        for (auto p = parent.getPolytope(d); p; ++p)
          p->getTransformation();
      }
      Point sample(*parent.getPolytope(d, 0), rc);
      expected = sample.getPhysicalCoordinates();
      for (size_t j = 0; j < d; ++j)
        for (size_t i = 0; i < parent.getSpaceDimension(); ++i)
          expectedJacobian.push_back(sample.getJacobian()(i, j));
      if (d == 0)
      {
        // The cell sharder requires positive-dimensional cells. Distribute
        // point entities directly through the same shard transport instead.
        for (int r = 0; r < world->size(); ++r)
        {
          Shard::Builder builder;
          builder.initialize(parent);
          builder.include({0, static_cast<Index>(r)}, Shard::State::Owned);
          sharder.getShards().push_back(builder.finalize());
        }
      }
      else
      {
        BalancedCompactPartitioner partitioner(parent);
        partitioner.partition(static_cast<size_t>(world->size()));
        sharder.shard(partitioner);
      }
    }
    // Broadcast before scatter waits for its sends: peers must be free to
    // enter gather and receive shards without relying on eager MPI buffering.
    boost::mpi::broadcast(*world, expected, 0);
    boost::mpi::broadcast(*world, expectedJacobian, 0);
    if (world->rank() == 0)
      sharder.scatter(0);
    auto mesh = sharder.gather(0);
    Location::AABB locator(mesh);
    locator.setProjectionPruning(true);
    size_t ownedHits = 0;
    for (auto p = mesh.getPolytope(d); p; ++p)
    {
      Point sample(*p, rc);
      auto hit = locator.locate(d, sample.getPhysicalCoordinates());
      if (!mesh.getShard().isOwned(d, p->getIndex()))
      {
        EXPECT_FALSE(hit);
        continue;
      }
      ASSERT_TRUE(hit) << "order=" << K << " embedded=" << embedded;
      ++ownedHits;
      EXPECT_EQ(&hit->getPolytope().getMesh(), &mesh);
      EXPECT_EQ(hit->getPolytope().getIndex(), p->getIndex());
      EXPECT_TRUE(mesh.isLocalPoint(*hit));
      EXPECT_LT((hit->getReferenceCoordinates() - rc).norm(), 1e-8);
      EXPECT_LT((hit->getJacobian() - sample.getJacobian()).norm(), 1e-8);
      if (mesh.getGlobalIndex(d, p->getIndex()) == 0)
      {
        EXPECT_LT((sample.getPhysicalCoordinates() - expected).norm(), 1e-12);
        for (size_t j = 0; j < d; ++j)
          for (size_t i = 0; i < mesh.getSpaceDimension(); ++i)
            EXPECT_NEAR(sample.getJacobian()(i, j),
              expectedJacobian[j * mesh.getSpaceDimension() + i], 1e-12);
      }
      // Field evaluation must accept the lifted point and use local entity IDs.
      Variational::P0<Real, MPIMesh> space(mesh);
      Variational::GridFunction field(space);
      field = Real(7);
      EXPECT_NEAR(field(*hit), Real(7), 1e-12);
    }
    EXPECT_EQ(boost::mpi::all_reduce(*world, ownedHits, std::plus<size_t>()),
      static_cast<size_t>(world->size()));
  }
}

TEST(MPI_Location_AABB, OwnedSubsetBeforeSearchAndFallback)
{
  Context::MPI context(*environment, *world);
  auto mesh = overlappingMesh(context);
  Location::AABB locator(mesh);
  locator.setExhaustiveFallback(true)
    .setProjectionPruning(true)
    .setTolerance(1e-9)
    .setReferenceTolerance(1e-9);
  Math::SpatialPoint x(2);
  x[0] = Real(0.25);
  x[1] = 0;
  auto hit = locator.locate(x);
  ASSERT_TRUE(hit);
  EXPECT_EQ(hit->getPolytope().getIndex(), 2);
  EXPECT_EQ(&hit->getPolytope().getMesh(), &mesh);
  EXPECT_EQ(mesh.getGlobalIndex(1, hit->getPolytope().getIndex()), 102);
  EXPECT_EQ(mesh.getLocalIndex(1, 102).value(), 2);
  EXPECT_FALSE(locator.locate(0, mesh.getShard().getVertexCoordinates(0)));
  x[1] = Real(0.01);
  EXPECT_FALSE(locator.locate(x));
}

TEST(MPI_Location_AABB, ConstructionAndQueriesOnOneRankOnly)
{
  Context::MPI context(*environment, *world);
  auto mesh = overlappingMesh(context);
  if (world->rank() == 0)
  {
    Location::AABB locator(mesh);
    Math::SpatialPoint x(2);
    x[0] = Real(0.5);
    x[1] = 0;
    EXPECT_TRUE(locator.locate(x));
    locator.setTolerance(1e-8).setExhaustiveFallback(true);
    EXPECT_TRUE(locator.locate(x));
  }
  world->barrier();
}

TEST(MPI_Location_AABB, EmptyRanksAndMPISubMesh)
{
  Context::MPI context(*environment, *world);
  auto parent = curvedMesh<3>(Polytope::Type::Triangle, true);
  Shard::Builder builder;
  builder.initialize(parent, 2);
  if (world->rank() == 0)
  {
    for (Index v = 0; v < 3; ++v)
      builder.include({0, v}, Shard::State::Owned);
    builder.include({2, 0}, Shard::State::Owned);
  }
  auto mesh = MPIMesh::Builder(context).initialize(builder.finalize()).finalize();
  SubMesh<Context::MPI>::Builder subBuilder;
  subBuilder.initialize(mesh);
  if (world->rank() == 0)
    subBuilder.include(2, 0);
  auto submesh = subBuilder.finalize();
  EXPECT_EQ(submesh.getDimension(), 2);
  EXPECT_EQ(submesh.getSpaceDimension(), parent.getSpaceDimension());
  Location::AABB locator(submesh);
  const auto rc = Polytope::Traits(Polytope::Type::Triangle).getCentroid();
  const auto x = Point(*parent.getCell(), rc).getPhysicalCoordinates();
  auto hit = locator.locate(x);
  ASSERT_EQ(hit.has_value(), world->rank() == 0);
  if (hit)
  {
    EXPECT_EQ(&hit->getPolytope().getMesh(), &submesh);
    EXPECT_LT((hit->getPhysicalCoordinates() - x).norm(), 1e-12);
  }
  Location::AABB meshLocator(mesh);
  EXPECT_EQ(meshLocator.locate(2, x).has_value(), world->rank() == 0);
}

TEST(MPI_Location_AABB, SharedBoundaryAndDimensionSpecificOwnership)
{
  Context::MPI context(*environment, *world);
  auto parent = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
  Shard::Builder builder;
  builder.initialize(parent);
  for (Index v = 0; v < 3; ++v)
  {
    const int owner = v == 2 && world->size() > 1 ? 1 : 0;
    const auto [i, inserted] = builder.include(
      {0, v}, world->rank() == owner ? Shard::State::Owned : Shard::State::Shared);
    if (world->rank() != owner)
      builder.setOwner(0, i, owner);
    else
      for (int r = 0; r < world->size(); ++r)
        if (r != owner)
          builder.halo(0, i, r);
  }
  for (Index c = 0; c < 2; ++c)
  {
    const int owner = static_cast<int>(c % world->size());
    const auto [i, inserted] = builder.include(
      {1, c}, world->rank() == owner ? Shard::State::Owned : Shard::State::Ghost);
    if (world->rank() != owner)
      builder.setOwner(1, i, owner);
    else
      for (int r = 0; r < world->size(); ++r)
        if (r != owner)
          builder.halo(1, i, r);
  }
  auto mesh = MPIMesh::Builder(context).initialize(builder.finalize()).finalize();
  Location::AABB locator(mesh);
  locator.setExhaustiveFallback(true);
  const auto& boundary = parent.getVertexCoordinates(1);
  const auto rc = Polytope::Traits(Polytope::Type::Segment).getCentroid();
  for (Index c = 0; c < 2; ++c)
  {
    const auto x = Point(*parent.getCell(c), rc).getPhysicalCoordinates();
    EXPECT_EQ(locator.locate(x).has_value(),
      world->rank() == static_cast<int>(c % world->size()));
  }
  EXPECT_EQ(locator.locate(1, boundary).has_value(), world->rank() < 2);
  EXPECT_EQ(locator.locate(0, boundary).has_value(), world->rank() == 0);
  for (int q = 0; q <= world->rank(); ++q)
    EXPECT_EQ(locator.locate(boundary).has_value(), world->rank() < 2);
  world->barrier();
}

class MPILocationGeometryTest : public ::testing::TestWithParam<Polytope::Type>
{};
TEST_P(MPILocationGeometryTest, TransfersAndLocatesAllOrders)
{
  for (bool embedded : {false, true})
  {
    checkDistributedGeometry<1>(GetParam(), embedded);
    checkDistributedGeometry<2>(GetParam(), embedded);
    checkDistributedGeometry<3>(GetParam(), embedded);
    checkDistributedGeometry<4>(GetParam(), embedded);
  }
}
TEST_P(MPILocationGeometryTest, TransfersCachedDefaultMaps)
{
  checkDistributedGeometry<1>(GetParam(), true, true);
}
INSTANTIATE_TEST_SUITE_P(
  AllGeometries, MPILocationGeometryTest, ::testing::ValuesIn(Polytope::Types));

int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator comm;
  environment = &env;
  world = &comm;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
