/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Unit tests for P0/P1/H1 spaces and native field restriction on
 * SubMesh<Context::MPI>.
 *
 * These tests verify that distributed P0 and P1 FES can be constructed on
 * distributed SubMeshes (both boundary skin and full-cell sub-regions)
 * and that the resulting DOF mappings are consistent across all MPI ranks:
 *
 *   - getSize() matches the expected global DOF count.
 *   - getVectorDimension() reports 1.
 *   - getOwnershipRange() width equals the number of locally owned DOFs.
 *   - getDOFs(D, i) and getGlobalIndex({D, i}, 0) are consistent.
 *   - All local polytopes have a valid (in-range) global DOF.
 *   - Owned DOFs are globally unique across all ranks.
 *
 * The separate H1 restriction matrix covers orders 1-6, real/complex scalar
 * vector and matrix fields, and all seven parent cell families, on affine and
 * exact quadratic domains. Full, boundary,
 * sparse and nested extractions retain exact logical ancestry. Every held
 * coefficient is compared with independent direct reference-polynomial interpolation;
 * physical point samples and a zero-field control provide additional oracles.
 * Distributed space setup and the global nonempty-coverage assertion are
 * collective. Native interpolation is also tested on rank zero alone;
 * protocol barriers are outside that operation.
 *
 * Legacy entries run at 1/2/3/4 ranks. The separate slow restriction entries
 * additionally run at 8 ranks to exercise empty holders.
 */

#include <set>
#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <type_traits>

#include <gtest/gtest.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/SubMesh.h>
#include <Rodin/MPI/Variational/P0g.h>
#include <Rodin/MPI/Variational/P0.h>
#include <Rodin/MPI/Variational/P1.h>
#include <Rodin/MPI/Variational/H1.h>
#include <Rodin/Variational.h>
#include "../../../../convergence/CurvedGeometry.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// ---------------------------------------------------------------------------
// Global MPI handles (initialized in main())
// ---------------------------------------------------------------------------
static boost::mpi::environment* g_env = nullptr;
static boost::mpi::communicator* g_world = nullptr;

namespace
{
  static Mesh<Context::Local> makeShardableMesh(
    Polytope::Type type, std::initializer_list<size_t> shape)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(type, shape);
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    if (D > 1)
    {
      mesh.getConnectivity().compute(D, 1);
      mesh.getConnectivity().compute(D - 1, 1);
    }
    return mesh;
  }

  static Mesh<Context::MPI> distributeFromRoot(const Context::MPI& ctx,
    Polytope::Type type, std::initializer_list<size_t> shape, bool curved = false)
  {
    const auto& comm = ctx.getCommunicator();
    Sharder<Context::MPI> sharder(ctx);
    if (comm.rank() == 0)
    {
      auto localMesh = makeShardableMesh(type, shape);
      if (curved)
      {
        Tests::Convergence::CurvedGeometry mapping(localMesh);
        mapping.template install<2>();
      }
      BalancedCompactPartitioner partitioner(localMesh);
      partitioner.partition(static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    return sharder.gather(0);
  }

  static const char* polytopeName(Polytope::Type type)
  {
    switch (type)
    {
      case Polytope::Type::Segment:
        return "Segment";
      case Polytope::Type::Triangle:
        return "Triangle";
      case Polytope::Type::Quadrilateral:
        return "Quadrilateral";
      case Polytope::Type::Tetrahedron:
        return "Tetrahedron";
      case Polytope::Type::Hexahedron:
        return "Hexahedron";
      case Polytope::Type::Pyramid:
        return "Pyramid";
      case Polytope::Type::Wedge:
        return "Wedge";
      default:
        return "Other";
    }
  }

  static SubMesh<Context::MPI> makeBoundarySubMesh(const Mesh<Context::MPI>& mesh)
  {
    const size_t faceDim = mesh.getDimension() - 1;
    SubMesh<Context::MPI>::Builder builder;
    builder.initialize(mesh);
    for (auto it = mesh.getBoundary(); it; ++it)
      builder.include(faceDim, it->getIndex());
    return builder.finalize();
  }

  static SubMesh<Context::MPI> makeCellSubMesh(const Mesh<Context::MPI>& mesh)
  {
    const size_t cellDim = mesh.getDimension();
    const auto& shard = mesh.getShard();
    SubMesh<Context::MPI>::Builder builder;
    builder.initialize(mesh);
    for (Index i = 0; i < static_cast<Index>(shard.getCellCount()); ++i)
    {
      if (shard.isOwned(cellDim, i))
        builder.include(cellDim, i);
    }
    return builder.finalize();
  }

  static SubMesh<Context::MPI> makeSparseBoundarySubMesh(
    const Mesh<Context::MPI>& mesh, const boost::mpi::communicator& comm)
  {
    const size_t faceDim = mesh.getDimension() - 1;
    const auto& shard = mesh.getShard();

    bool hasOwnedBoundary = false;
    Index selectedFace = std::numeric_limits<Index>::max();
    for (auto it = mesh.getBoundary(); it; ++it)
    {
      const Index faceIdx = it->getIndex();
      if (shard.isOwned(faceDim, faceIdx))
      {
        hasOwnedBoundary = true;
        selectedFace = faceIdx;
        break;
      }
    }

    std::vector<int> hasBoundaryByRank;
    boost::mpi::all_gather(comm, hasOwnedBoundary ? 1 : 0, hasBoundaryByRank);

    int selectedRank = -1;
    for (int r = 0; r < static_cast<int>(hasBoundaryByRank.size()); ++r)
    {
      if (hasBoundaryByRank[static_cast<size_t>(r)])
      {
        selectedRank = r;
        break;
      }
    }
    assert(selectedRank >= 0);

    SubMesh<Context::MPI>::Builder builder;
    builder.initialize(mesh);
    if (comm.rank() == selectedRank)
      builder.include(faceDim, selectedFace);
    return builder.finalize();
  }
}

namespace Rodin::Tests::Unit
{
  /** @brief Reference polynomial restriction with exact logical provenance.
   * H1 orders one through six and all six scalar/vector/matrix value ranges are
   * checked on full, boundary, sparse and nested distributed SubMeshes.
   * The floating-point budget concerns field values and DOF-functional
   * evaluation, never entity correspondence. Space construction is collective; native interpolation
   * is deliberately also exercised on rank zero alone.
   * On curved domains an independent analytic inverse evaluates the reference
   * polynomial; it is never used to identify entities or DOFs. The quadratic
   * map is installed before partitioning and must survive child extraction.
   * With ambient dimension @f$d@f$ and reference coordinates @f$x@f$,
   * the degree-@f$K@f$ mixed terms use
   * @f$s_r=\sum_j(j+1)x_j/[d(d+1)/2]@f$ and
   * @f$s_i=\sum_j(j+2)x_j/[d(d+3)/2]@f$.
   * For @f$K\geq2@f$, these enter as @f$s_r^K,s_i^K@f$;
   * @f$K=1@f$ uses affine data. Their unit-box magnitudes stay bounded
   * independently of @f$K@f$ while exciting mixed high-order modes.
   */
  TEST(MPISubMeshH1, ExplicitRestriction_AllOrdersAllRangesAllGeometries)
  {
    const auto& world = *g_world;
    Context::MPI context(*g_env, world);
    constexpr Real PolynomialTolerance = 1e-10;
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(polytopeName(geometry));
      for (bool curved : {false, true})
      {
        SCOPED_TRACE(::testing::Message() << "curved=" << curved);
        const size_t dimension = Polytope::Traits(geometry).getDimension();
        auto parent = dimension == 1 ? distributeFromRoot(context, geometry, {2}, curved)
          : dimension == 2 ? distributeFromRoot(context, geometry, {2, 2}, curved)
                           : distributeFromRoot(context, geometry, {2, 2, 2}, curved);
        const auto check = [&](const SubMesh<Context::MPI>& sub) {
          const Mesh<Context::MPI>& mesh = sub;
          size_t localOwned = 0;
          for (auto cell = mesh.getCell(); cell; ++cell)
          {
            if (mesh.getShard().isOwned(mesh.getDimension(), cell->getIndex()))
              ++localOwned;
          }
        // This is a global coverage assertion, not a local evaluation query.
          EXPECT_GT(boost::mpi::all_reduce(world, localOwned, std::plus<size_t>()), 0u);
          const auto range = [&]<size_t K, class Range>() {
            using Space = H1<K, Range, Mesh<Context::MPI>>;
            using Scalar = typename Space::ScalarType;
            constexpr bool vector = std::is_same_v<Range, Math::SpatialVector<Scalar>>;
            constexpr bool matrix = std::is_same_v<Range, Math::SpatialMatrix<Scalar>>;
            SCOPED_TRACE((::testing::Message()
              << "degree=" << K << " vector=" << vector << " matrix=" << matrix
              << " complex=" << std::is_same_v<Scalar, Complex>));
            const auto makeSpace = [](const Mesh<Context::MPI>& mesh) {
              if constexpr (matrix)
                return Space(std::integral_constant<size_t, K>{}, mesh, 2, 3);
              else if constexpr (vector)
                return Space(std::integral_constant<size_t, K>{}, mesh, 3);
              else
                return Space(std::integral_constant<size_t, K>{}, mesh);
            };
            const auto sourceSpace = makeSpace(parent), targetSpace = makeSpace(mesh);
            GridFunction source(sourceSpace), restricted(targetSpace),
              oracle(targetSpace), rejected(targetSpace);
            rejected.getData().setZero();
            const auto scalar = [curved](const Point& point, size_t component) {
              Math::SpatialPoint x = point.getPhysicalCoordinates();
              if (curved)
              {
                // Independent analytic inverse defines the field, not correspondence.
                if (x.size() == 1)
                  x(0) = 2 * x(0) / (1 + std::sqrt(1 + Real(0.4) * x(0)));
                else
                  x(x.size() - 1) -= Real(0.1) * x(0) * x(0);
              }
              Real re = 1, im = 2;
              for (Index j = 0; j < x.size(); ++j)
              {
                re += Real(j + 1) * x(j);
                im += Real(j + 2) * x(j);
              }
              if constexpr (K >= 2)
              {
                // Normalized mixed degree-K polynomials excite higher modes
                // without growing the dimensionless field scale with K.
                const Real dim = Real(x.size());
                re += std::pow((re - 1) / (dim * (dim + 1) / 2), Real(K));
                im += std::pow((im - 2) / (dim * (dim + 3) / 2), Real(K));
              }
              if constexpr (std::is_same_v<Scalar, Complex>)
                return Real(component + 1) * Complex(re, im);
              else
                return Real(component + 1) * re;
            };
            const auto data = [&] {
              if constexpr (matrix)
                return MatrixFunction(size_t{2}, size_t{3}, [scalar](const Point& point) {
                  Math::SpatialMatrix<Scalar> value(2, 3);
                  for (size_t row = 0; row < 2; ++row)
                  {
                    for (size_t column = 0; column < 3; ++column)
                      value(row, column) = scalar(point, 3 * row + column);
                  }
                  return value;
                });
              else if constexpr (vector)
                return VectorFunction(size_t{3}, [scalar](const Point& point) {
                  Math::SpatialVector<Scalar> value(3);
                  for (size_t component = 0; component < 3; ++component)
                    value(component) = scalar(point, component);
                  return value;
                });
              else if constexpr (std::is_same_v<Scalar, Complex>)
                return ComplexFunction(
                  [scalar](const Point& point) { return scalar(point, 0); });
              else
                return RealFunction(
                  [scalar](const Point& point) { return scalar(point, 0); });
            }();
            source.project(data);
            oracle.project(data);
            if (world.rank() == 0)
            {
              restricted.project(source);
              // Check the one-rank result before the subsequent all-rank call
              // can overwrite it. Metadata queries and held-DOF access are local.
              EXPECT_EQ(&restricted.getFiniteElementSpace(), &targetSpace);
              EXPECT_EQ(restricted.getData().size(), targetSpace.getSize());
              for (auto cell = mesh.getCell(); cell; ++cell)
              {
                const IndexArray dofs =
                  targetSpace.getDOFs(mesh.getDimension(), cell->getIndex());
                for (Index local = 0; local < static_cast<size_t>(dofs.size()); ++local)
                {
                  EXPECT_LT(std::abs(restricted[dofs(local)] - oracle[dofs(local)]),
                    PolynomialTolerance);
                }
                // Point evaluation must also remain local while peers wait.
                const Point point(
                  *cell, Polytope::Traits(cell->getGeometry()).getCentroid());
                const auto actual = std::as_const(restricted)(point),
                           expected = data(point);
                if constexpr (matrix)
                {
                  EXPECT_EQ(actual.rows(), 2);
                  EXPECT_EQ(actual.cols(), 3);
                  for (size_t row = 0; row < 2; ++row)
                  {
                    for (size_t column = 0; column < 3; ++column)
                    {
                      EXPECT_LT(std::abs(actual(row, column) - expected(row, column)),
                        PolynomialTolerance);
                    }
                  }
                }
                else if constexpr (vector)
                {
                  for (size_t component = 0; component < 3; ++component)
                  {
                    EXPECT_LT(std::abs(actual(component) - expected(component)),
                      PolynomialTolerance);
                  }
                }
                else
                  EXPECT_LT(std::abs(actual - expected), PolynomialTolerance);
              }
            }
            // Test-protocol synchronization, outside the local value operation.
            world.barrier();
            restricted.project(source);
            EXPECT_EQ(&restricted.getFiniteElementSpace(), &targetSpace);
            EXPECT_EQ(restricted.getData().size(), targetSpace.getSize());
            const size_t d = mesh.getDimension();
            for (auto cell = mesh.getCell(); cell; ++cell)
            {
              const auto& ancestry = sub.getPolytopeMap(d).left;
              ASSERT_LT(cell->getIndex(), ancestry.size());
              Index ancestor = ancestry[cell->getIndex()];
              const MeshBase* immediate = &sub.getParent();
              while (immediate->isSubMesh())
              {
                const auto& upper = immediate->asSubMesh();
                const auto& map = upper.getPolytopeMap(d).left;
                ASSERT_LT(ancestor, map.size());
                ancestor = map[ancestor];
                immediate = &upper.getParent();
              }
              ASSERT_EQ(immediate, &parent);
              const auto original = parent.getPolytope(d, ancestor);
              ASSERT_TRUE(original);
              EXPECT_EQ(original->getGeometry(), cell->getGeometry());
              EXPECT_EQ(cell->getTransformation().getFactorOrder(),
                original->getTransformation().getFactorOrder());
              if (curved && d > 0)
              {
                EXPECT_EQ(cell->getTransformation().getFactorOrder(), 2u);
              }
              const IndexArray dofs = targetSpace.getDOFs(d, cell->getIndex());
              for (Index local = 0; local < static_cast<size_t>(dofs.size()); ++local)
              {
                EXPECT_EQ(
                  dofs(local), targetSpace.getGlobalIndex({d, cell->getIndex()}, local));
                // Every held functional is checked: sparse point samples alone
                // cannot certify a high-order polynomial restriction.
                EXPECT_LT(std::abs(restricted[dofs(local)] - oracle[dofs(local)]),
                  PolynomialTolerance);
              }
              const Polytope::Traits traits(cell->getGeometry());
              const auto sample = [&](const Math::SpatialPoint& coordinates) {
                const Point point(*cell, coordinates);
                const auto actual = restricted(point), expected = data(point),
                           wrong = rejected(point);
                if constexpr (matrix)
                {
                  EXPECT_EQ(actual.rows(), 2);
                  EXPECT_EQ(actual.cols(), 3);
                  for (size_t row = 0; row < 2; ++row)
                  {
                    for (size_t column = 0; column < 3; ++column)
                    {
                      EXPECT_LT(std::abs(actual(row, column) - expected(row, column)),
                        PolynomialTolerance);
                      EXPECT_GT(std::abs(wrong(row, column) - expected(row, column)),
                        PolynomialTolerance);
                    }
                  }
                }
                else if constexpr (vector)
                {
                  for (size_t component = 0; component < 3; ++component)
                  {
                    EXPECT_LT(std::abs(actual(component) - expected(component)),
                      PolynomialTolerance);
                    EXPECT_GT(std::abs(wrong(component) - expected(component)),
                      PolynomialTolerance);
                  }
                }
                else
                {
                  EXPECT_LT(std::abs(actual - expected), PolynomialTolerance);
                  EXPECT_GT(std::abs(wrong - expected), PolynomialTolerance);
                }
              };
              sample(traits.getCentroid());
              for (size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
              {
                sample(Math::SpatialPoint(
                  (traits.getCentroid() + traits.getVertex(vertex)) / 2));
              }
            }
          };
          const auto order = [&]<size_t K>() {
            range.template operator()<K, Real>();
            range.template operator()<K, Complex>();
            range.template operator()<K, Math::SpatialVector<Real>>();
            range.template operator()<K, Math::SpatialVector<Complex>>();
            range.template operator()<K, Math::SpatialMatrix<Real>>();
            range.template operator()<K, Math::SpatialMatrix<Complex>>();
          };
          order.template operator()<1>();
          order.template operator()<2>();
          order.template operator()<3>();
          order.template operator()<4>();
          order.template operator()<5>();
          order.template operator()<6>();
        };
        const auto cells = makeCellSubMesh(parent);
        check(cells);
        const auto skin = makeBoundarySubMesh(parent);
        check(skin);
        SubMesh<Context::MPI>::Builder builder;
        builder.initialize(parent);
        for (auto cell = parent.getCell(); cell; ++cell)
        {
          if (parent.getShard().isOwned(dimension, cell->getIndex()) &&
            parent.getGlobalIndex(dimension, cell->getIndex()) % 2 == 0)
            builder.include(dimension, cell->getIndex());
        }
        const auto sparse = builder.finalize();
        check(sparse);
        const auto nested = makeCellSubMesh(sparse);
        check(nested);
      }
    }
  }

  /// @brief Verifies empty ranks construct P0 P1 H1 triangle boundary for MPI sub mesh sparse selection by checking exact expected values, MPI behavior.
  TEST(MPISubMeshSparseSelection, EmptyRanks_ConstructP0P1H1_TriangleBoundary)
  {
    const auto& world = *g_world;
    if (world.size() < 2)
      GTEST_SKIP() << "Requires at least two MPI ranks to exercise empty submesh ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeSparseBoundarySubMesh(mesh, world);

    const size_t faceDim = mesh.getDimension() - 1;
    ASSERT_EQ(sub.getDimension(), faceDim);
    ASSERT_EQ(sub.getPolytopeCount(faceDim), 1u);

    P0<Real, Mesh<Context::MPI>> p0(sub);
    EXPECT_EQ(p0.getSize(), 1u);

    P1<Real, Mesh<Context::MPI>> p1(sub);
    const size_t globalVertices = sub.getPolytopeCount(0);
    EXPECT_EQ(p1.getSize(), globalVertices);

    H1<2, Real, Mesh<Context::MPI>> h1(std::integral_constant<size_t, 2>{}, sub);
    EXPECT_EQ(h1.getSize(), globalVertices + 1u);

    Index begin = 0;
    Index end = 0;
    h1.getOwnershipRange(begin, end);
    const size_t localOwned = static_cast<size_t>(end - begin);
    size_t globalOwned = 0;
    boost::mpi::all_reduce(world, localOwned, globalOwned, std::plus<size_t>());
    EXPECT_EQ(globalOwned, h1.getSize());
  }

  /// @brief Verifies boundary sub mesh scalar global constant DO fs triangle for MPI sub mesh P0 g by checking exact expected values, MPI behavior.
  TEST(MPISubMeshP0g, BoundarySubMesh_ScalarGlobalConstantDOFs_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeBoundarySubMesh(mesh);

    P0g<Real, Mesh<Context::MPI>> fes(sub);
    EXPECT_EQ(fes.getSize(), 1u);
    EXPECT_EQ(fes.getVectorDimension(), 1u);

    Index begin = 0;
    Index end = 0;
    fes.getOwnershipRange(begin, end);
    if (world.rank() == 0)
    {
      EXPECT_EQ(begin, Index(0));
      EXPECT_EQ(end, Index(1));
    }
    else
    {
      EXPECT_EQ(begin, Index(1));
      EXPECT_EQ(end, Index(1));
    }

    const size_t d = sub.getDimension();
    const auto& shard = sub.getShard();
    for (Index i = 0; i < static_cast<Index>(shard.getPolytopeCount(d)); ++i)
    {
      const auto& dofs = fes.getDOFs(d, i);
      ASSERT_EQ(dofs.size(), 1u);
      EXPECT_EQ(dofs[0], Index(0));
      EXPECT_EQ(fes.getGlobalIndex({d, i}, 0), Index(0));
    }
  }

  /// @brief Verifies cell sub mesh vector global constant DO fs triangle for MPI sub mesh P0 g by checking exact expected values, MPI behavior.
  TEST(MPISubMeshP0g, CellSubMesh_VectorGlobalConstantDOFs_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeCellSubMesh(mesh);

    const size_t vdim = 3;
    P0g<Math::SpatialVector<Real>, Mesh<Context::MPI>> fes(sub, vdim);
    EXPECT_EQ(fes.getSize(), vdim);
    EXPECT_EQ(fes.getVectorDimension(), vdim);

    Index begin = 0;
    Index end = 0;
    fes.getOwnershipRange(begin, end);
    if (world.rank() == 0)
    {
      EXPECT_EQ(begin, Index(0));
      EXPECT_EQ(end, static_cast<Index>(vdim));
    }
    else
    {
      EXPECT_EQ(begin, static_cast<Index>(vdim));
      EXPECT_EQ(end, static_cast<Index>(vdim));
    }

    const size_t d = sub.getDimension();
    const auto& shard = sub.getShard();
    for (Index i = 0; i < static_cast<Index>(shard.getPolytopeCount(d)); ++i)
    {
      const auto& dofs = fes.getDOFs(d, i);
      ASSERT_EQ(dofs.size(), vdim);
      for (size_t k = 0; k < vdim; ++k)
      {
        EXPECT_EQ(dofs[k], static_cast<Index>(k));
        EXPECT_EQ(
          fes.getGlobalIndex({d, i}, static_cast<Index>(k)), static_cast<Index>(k));
      }
    }
  }

  /// @brief Verifies sparse boundary sub mesh constructs on empty ranks for MPI sub mesh P0 g by checking exact expected values, MPI behavior.
  TEST(MPISubMeshP0g, SparseBoundarySubMesh_ConstructsOnEmptyRanks)
  {
    const auto& world = *g_world;
    if (world.size() < 2)
      GTEST_SKIP() << "Requires at least two MPI ranks to exercise empty submesh ranks.";
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeSparseBoundarySubMesh(mesh, world);

    const size_t d = sub.getDimension();
    ASSERT_EQ(sub.getPolytopeCount(d), 1u);

    P0g<Real, Mesh<Context::MPI>> fes(sub);
    EXPECT_EQ(fes.getSize(), 1u);

    const auto& shard = sub.getShard();
    for (Index i = 0; i < static_cast<Index>(shard.getPolytopeCount(d)); ++i)
    {
      const auto& dofs = fes.getDOFs(d, i);
      ASSERT_EQ(dofs.size(), 1u);
      EXPECT_EQ(dofs[0], Index(0));
    }
  }

  /// @brief Verifies boundary sub mesh quadratic H1 all 3 D for MPI sub mesh H1 by checking exact expected values, MPI behavior.
  TEST(MPISubMeshH1, BoundarySubMesh_QuadraticH1_All3D)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    for (auto type : {Polytope::Type::Tetrahedron, Polytope::Type::Hexahedron,
           Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(polytopeName(type));
      auto mesh = distributeFromRoot(ctx, type, {4, 3, 3});
      auto sub = makeBoundarySubMesh(mesh);

      ASSERT_EQ(sub.getDimension(), 2u);

      H1<2, Real, Mesh<Context::MPI>> h1(std::integral_constant<size_t, 2>{}, sub);
      EXPECT_GT(h1.getSize(), 0u);

      const auto& shard = sub.getShard();
      const auto& connectivity = sub.getConnectivity();
      for (Index i = 0; i < static_cast<Index>(shard.getPolytopeCount(2)); ++i)
      {
        const auto& dofs = h1.getDOFs(2, i);
        const auto geometry = connectivity.getGeometry(2, i);
        EXPECT_EQ(dofs.size(), RealH1Element<2>(geometry).getCount());
        for (const Index dof : dofs)
          EXPECT_LT(dof, static_cast<Index>(h1.getSize()));
      }
    }
  }

  // ==========================================================================
  // Group 1 — P0 FES on boundary SubMesh<Context::MPI>
  // ==========================================================================

  /**
   * @brief Scalar P0 FES on boundary SubMesh: global DOF count equals the
   *        number of boundary faces.
   */
  TEST(MPIP0FESSubMesh, BoundarySubMesh_GetSize_EqualsFaceCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);

    // P0 DOFs = number of faces (cells of the SubMesh).
    const size_t globalFaces = sub.getPolytopeCount(sub.getDimension());
    EXPECT_EQ(fes.getSize(), globalFaces);
  }

  /**
   * @brief Scalar P0 on boundary SubMesh has vector dimension 1.
   */
  TEST(MPIP0FESSubMesh, BoundarySubMesh_VectorDimension_IsOne_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);
    EXPECT_EQ(fes.getVectorDimension(), 1u);
  }

  /**
   * @brief Ownership range of P0 on boundary SubMesh has width equal to
   *        the number of locally owned boundary faces.
   */
  TEST(MPIP0FESSubMesh, BoundarySubMesh_OwnershipRange_MatchesOwnedFaces_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard = sub.getShard();
    const size_t fDim = sub.getDimension();
    const size_t nFace = shard.getCellCount();

    size_t ownedCount = 0;
    for (size_t i = 0; i < nFace; ++i)
    {
      if (shard.isOwned(fDim, i))
        ++ownedCount;
    }

    Index begin = 0, end = 0;
    fes.getOwnershipRange(begin, end);
    EXPECT_EQ(static_cast<size_t>(end - begin), ownedCount);
  }

  /**
   * @brief getDOFs and getGlobalIndex are consistent for P0 on boundary SubMesh.
   *
   * For P0 each element has exactly one DOF.  getDOFs(D, i)[0] must equal
   * getGlobalIndex({D, i}, 0).
   */
  TEST(MPIP0FESSubMesh, BoundarySubMesh_GetDOFs_ConsistentWithGetGlobalIndex_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard = sub.getShard();
    const size_t fDim = sub.getDimension();
    const size_t nFace = shard.getCellCount();

    for (size_t i = 0; i < nFace; ++i)
    {
      const auto& dofs = fes.getDOFs(fDim, static_cast<Index>(i));
      ASSERT_EQ(dofs.size(), 1u);
      const Index via_dofs = dofs[0];
      const Index via_global = fes.getGlobalIndex({fDim, static_cast<Index>(i)}, 0);
      EXPECT_EQ(via_dofs, via_global) << "Mismatch for face " << i;
    }
  }

  /**
   * @brief All local boundary faces have a valid global P0 DOF index.
   */
  TEST(MPIP0FESSubMesh, BoundarySubMesh_AllLocalFaces_HaveValidGlobalDOF_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeBoundarySubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard = sub.getShard();
    const size_t fDim = sub.getDimension();
    const size_t nFace = shard.getCellCount();

    // The global size is established during collective space construction;
    // getSize() only reads cached metadata and is safe on an individual rank.
    const size_t globalSize = fes.getSize();

    for (size_t i = 0; i < nFace; ++i)
    {
      const auto& dofs = fes.getDOFs(fDim, static_cast<Index>(i));
      for (const Index d : dofs)
      {
        EXPECT_LT(d, globalSize)
          << "Global P0 DOF " << d << " out of range for local face " << i;
      }
    }
  }

  // ==========================================================================
  // Group 2 — P0 FES on full-cell SubMesh<Context::MPI>
  // ==========================================================================

  /**
   * @brief P0 FES on full-cell SubMesh: global DOF count equals the parent
   *        mesh global cell count (SubMesh covers the whole domain).
   */
  TEST(MPIP0FESSubMesh, CellSubMesh_GetSize_EqualsParentCellCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeCellSubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);

    // Cell SubMesh is built from all owned cells, so global size = global cell count.
    // Collective call — issued by all ranks.
    const size_t parentCells = mesh.getCellCount();
    EXPECT_EQ(fes.getSize(), parentCells);
  }

  /**
   * @brief Owned P0 DOFs on the full-cell SubMesh are globally unique.
   *
   * Gather all owned DOFs on rank 0 and verify no duplicates appear.
   */
  TEST(MPIP0FESSubMesh, CellSubMesh_OwnedDOFs_GloballyUnique_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeCellSubMesh(mesh);

    P0<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard = sub.getShard();
    const size_t cellDim = sub.getDimension();
    const size_t nCells = shard.getCellCount();

    std::vector<Index> ownedDofs;
    for (size_t i = 0; i < nCells; ++i)
    {
      if (shard.isOwned(cellDim, i))
      {
        const auto& dofs = fes.getDOFs(cellDim, static_cast<Index>(i));
        ASSERT_EQ(dofs.size(), 1u);
        ownedDofs.push_back(dofs[0]);
      }
    }

    // Collective call — must be issued by all ranks.
    const size_t globalCells = mesh.getCellCount();

    std::vector<std::vector<Index>> allDofs;
    boost::mpi::gather(world, ownedDofs, allDofs, 0);

    if (world.rank() == 0)
    {
      std::vector<Index> combined;
      for (const auto& v : allDofs)
        combined.insert(combined.end(), v.begin(), v.end());

      std::sort(combined.begin(), combined.end());
      const size_t uniqueCount = static_cast<size_t>(
        std::unique(combined.begin(), combined.end()) - combined.begin());

      EXPECT_EQ(uniqueCount, combined.size())
        << "Duplicate global P0 DOF indices on cell SubMesh.";
      EXPECT_EQ(combined.size(), globalCells);
    }
  }

  // ==========================================================================
  // Group 3 — 3D mesh variants
  // ==========================================================================

  /**
   * @brief P0 on boundary SubMesh of every 3D cell mesh: size equals face count.
   */
  TEST(MPIP0FESSubMesh, BoundarySubMesh_GetSize_EqualsFaceCount_All3D)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    for (auto type : {Polytope::Type::Tetrahedron, Polytope::Type::Hexahedron,
           Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(polytopeName(type));
      auto mesh = distributeFromRoot(ctx, type, {4, 3, 3});
      auto sub = makeBoundarySubMesh(mesh);

      P0<Real, Mesh<Context::MPI>> fes(sub);

      const size_t globalFaces = sub.getPolytopeCount(sub.getDimension());
      EXPECT_EQ(fes.getSize(), globalFaces);
    }
  }

  // ==========================================================================
  // Group 4 — P1 FES on boundary SubMesh<Context::MPI>
  // ==========================================================================

  /**
   * @brief Scalar P1 FES on boundary SubMesh: global DOF count equals the
   *        number of boundary vertices.
   */
  TEST(MPIP1FESSubMesh, BoundarySubMesh_GetSize_EqualsVertexCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const size_t globalVerts = sub.getPolytopeCount(0);
    EXPECT_EQ(fes.getSize(), globalVerts);
  }

  /**
   * @brief Ownership range of P1 on boundary SubMesh has width equal to
   *        the number of locally owned boundary vertices.
   */
  TEST(MPIP1FESSubMesh, BoundarySubMesh_OwnershipRange_MatchesOwnedVertices_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard = sub.getShard();
    const size_t nVert = shard.getVertexCount();

    size_t ownedCount = 0;
    for (size_t i = 0; i < nVert; ++i)
    {
      if (shard.isOwned(0, i))
        ++ownedCount;
    }

    Index begin = 0, end = 0;
    fes.getOwnershipRange(begin, end);
    EXPECT_EQ(static_cast<size_t>(end - begin), ownedCount);
  }

  /**
   * @brief All local boundary vertices have a valid global P1 DOF index.
   */
  TEST(MPIP1FESSubMesh, BoundarySubMesh_AllLocalVertices_HaveValidGlobalDOF_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard = sub.getShard();
    const size_t nVert = shard.getVertexCount();

    const size_t globalSize = fes.getSize();

    for (size_t i = 0; i < nVert; ++i)
    {
      const auto& dofs = fes.getDOFs(0, static_cast<Index>(i));
      for (const Index d : dofs)
      {
        EXPECT_LT(d, globalSize)
          << "Global P1 DOF " << d << " out of range for local vertex " << i;
      }
    }
  }

  // ==========================================================================
  // Group 5 — P1 FES on full-cell SubMesh<Context::MPI>
  // ==========================================================================

  /**
   * @brief P1 FES on full-cell SubMesh: global DOF count equals the parent
   *        mesh global vertex count.
   */
  TEST(MPIP1FESSubMesh, CellSubMesh_GetSize_EqualsParentVertexCount_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeCellSubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const size_t parentVerts = mesh.getVertexCount();
    EXPECT_EQ(fes.getSize(), parentVerts);
  }

  /**
   * @brief Owned P1 DOFs on the full-cell SubMesh are globally unique.
   *
   * Gather all owned DOFs on rank 0 and verify no duplicates appear.
   */
  TEST(MPIP1FESSubMesh, CellSubMesh_OwnedDOFs_GloballyUnique_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeCellSubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard = sub.getShard();
    const size_t nVerts = shard.getVertexCount();

    std::vector<Index> ownedDofs;
    for (size_t i = 0; i < nVerts; ++i)
    {
      if (shard.isOwned(0, i))
      {
        const auto& dofs = fes.getDOFs(0, static_cast<Index>(i));
        for (const Index d : dofs)
          ownedDofs.push_back(d);
      }
    }

    const size_t globalVerts = mesh.getVertexCount();

    std::vector<std::vector<Index>> allDofs;
    boost::mpi::gather(world, ownedDofs, allDofs, 0);

    if (world.rank() == 0)
    {
      std::vector<Index> combined;
      for (const auto& v : allDofs)
        combined.insert(combined.end(), v.begin(), v.end());

      std::sort(combined.begin(), combined.end());
      const size_t uniqueCount = static_cast<size_t>(
        std::unique(combined.begin(), combined.end()) - combined.begin());

      EXPECT_EQ(uniqueCount, combined.size())
        << "Duplicate global P1 DOF indices on cell SubMesh.";
      EXPECT_EQ(combined.size(), globalVerts);
    }
  }

  // ==========================================================================
  // Group 6 — P1 FES on boundary SubMesh: owned DOF uniqueness
  // ==========================================================================

  /**
   * @brief Owned P1 DOFs on boundary SubMesh are globally unique and cover
   *        the expected range exactly.
   *
   * This is a regression test for the MPI tag-collision bug: before the fix,
   * P0 cell SubMesh sends with tag=0 would accumulate and be consumed by a
   * subsequent P1 boundary SubMesh irecv on the same tag, leaving some
   * Shared vertices with invalid DOFs.
   */
  TEST(MPIP1FESSubMesh, BoundarySubMesh_OwnedDOFs_GloballyUnique_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});
    auto sub = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard = sub.getShard();
    const size_t nVerts = shard.getVertexCount();

    std::vector<Index> ownedDofs;
    for (size_t i = 0; i < nVerts; ++i)
    {
      if (shard.isOwned(0, i))
      {
        const auto& dofs = fes.getDOFs(0, static_cast<Index>(i));
        for (const Index d : dofs)
          ownedDofs.push_back(d);
      }
    }

    // Collective gather — must be called by all ranks.
    const size_t globalBoundaryVerts = sub.getVertexCount();

    std::vector<std::vector<Index>> allDofs;
    boost::mpi::gather(world, ownedDofs, allDofs, 0);

    if (world.rank() == 0)
    {
      std::vector<Index> combined;
      for (const auto& v : allDofs)
        combined.insert(combined.end(), v.begin(), v.end());

      std::sort(combined.begin(), combined.end());
      const size_t uniqueCount = static_cast<size_t>(
        std::unique(combined.begin(), combined.end()) - combined.begin());

      EXPECT_EQ(uniqueCount, combined.size())
        << "Duplicate global P1 DOF indices on boundary SubMesh.";
      EXPECT_EQ(combined.size(), globalBoundaryVerts)
        << "Owned P1 DOF count does not match global boundary vertex count.";
    }
  }

  // ==========================================================================
  // Group 7 — Regression: P0 cell SubMesh then P1 boundary SubMesh (tag-leakage)
  // ==========================================================================

  /**
   * @brief Regression test for MPI message tag collision.
   *
   * Previously, P0 on a cell SubMesh used tag=0 for its DOF exchange and
   * could leave unmatched sends in the MPI buffer.  A subsequent P1 boundary
   * SubMesh construction would irecv on tag=0 and consume those stale P0
   * messages, corrupting Shared-vertex DOF assignments.
   *
   * After the fix (symmetric drain, distinct tags), the P1 boundary SubMesh
   * must produce valid DOFs even when constructed after P0 cell and P0
   * boundary SubMesh spaces.
   */
  TEST(
    MPIP1FESSubMesh, Regression_P0CellThenP1Boundary_AllLocalVertices_HaveValidGlobalDOF)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {4, 4});

    // Build both SubMeshes.
    auto cellSub = makeCellSubMesh(mesh);
    auto boundarySub = makeBoundarySubMesh(mesh);

    // Construct P0 on cell SubMesh first (this was the source of stale sends).
    P0<Real, Mesh<Context::MPI>> p0cell(cellSub);
    (void)p0cell;

    // Construct P0 on boundary SubMesh (additional potential source of stale sends).
    P0<Real, Mesh<Context::MPI>> p0bnd(boundarySub);
    (void)p0bnd;

    // Now construct P1 on boundary SubMesh.  After the symmetric-drain fix,
    // this must not consume any stale P0 messages and must produce valid DOFs.
    P1<Real, Mesh<Context::MPI>> p1bnd(boundarySub);

    const auto& shard = boundarySub.getShard();
    const size_t nVerts = shard.getVertexCount();
    const size_t globalSize = p1bnd.getSize();

    for (size_t i = 0; i < nVerts; ++i)
    {
      const auto& dofs = p1bnd.getDOFs(0, static_cast<Index>(i));
      for (const Index d : dofs)
      {
        EXPECT_LT(d, globalSize)
          << "P1 boundary SubMesh DOF " << d << " out of range for local vertex " << i
          << " (constructed after P0 cell and boundary SubMesh spaces).";
      }
    }
  }

  /**
   * @brief Same regression test over all supported 3D cell meshes.
   */
  TEST(MPIP1FESSubMesh,
    Regression_P0CellThenP1Boundary_AllLocalVertices_HaveValidGlobalDOF_All3D)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);
    for (auto type : {Polytope::Type::Tetrahedron, Polytope::Type::Hexahedron,
           Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(polytopeName(type));
      auto mesh = distributeFromRoot(ctx, type, {4, 3, 3});

      auto cellSub = makeCellSubMesh(mesh);
      auto boundarySub = makeBoundarySubMesh(mesh);

      P0<Real, Mesh<Context::MPI>> p0cell(cellSub);
      (void)p0cell;

      P0<Real, Mesh<Context::MPI>> p0bnd(boundarySub);
      (void)p0bnd;

      P1<Real, Mesh<Context::MPI>> p1bnd(boundarySub);

      const auto& shard = boundarySub.getShard();
      const size_t nVerts = shard.getVertexCount();
      const size_t globalSize = p1bnd.getSize();

      for (size_t i = 0; i < nVerts; ++i)
      {
        const auto& dofs = p1bnd.getDOFs(0, static_cast<Index>(i));
        for (const Index d : dofs)
        {
          EXPECT_LT(d, globalSize)
            << "P1 boundary SubMesh DOF " << d << " out of range for local vertex " << i;
        }
      }
    }
  }

  // ==========================================================================
  // Group 8 — Regression: SubMesh orphaned vertex (2-round ownership fix)
  // ==========================================================================

  /**
   * @brief Regression test for the SubMesh orphaned-vertex ownership bug.
   *
   * A vertex v with Shared(owner=A) in the volume partition is added to the
   * boundary SubMesh on ranks B (and possibly C,D,...) when those ranks have
   * boundary faces adjacent to v.  If rank A has no boundary face adjacent
   * to v, the SubMesh finalize() 2-round exchange must promote the
   * minimum-rank querier to Owned so that P1 can assign a valid global DOF.
   *
   * Passing this test requires the 2-round ownership fix in SubMesh::finalize.
   */
  TEST(
    MPIP1FESSubMesh, Regression_OrphanVertex_AllLocalVertices_HaveValidGlobalDOF_Triangle)
  {
    const auto& world = *g_world;
    if (world.size() > 4)
      GTEST_SKIP() << "Test designed for at most 4 MPI ranks.";

    Context::MPI ctx(*g_env, world);

    // A finer 8×8 mesh increases the probability of hitting corner vertices
    // where the volume-partition owner is not adjacent to any boundary face.
    auto mesh = distributeFromRoot(ctx, Polytope::Type::Triangle, {8, 8});
    auto sub = makeBoundarySubMesh(mesh);

    P1<Real, Mesh<Context::MPI>> fes(sub);

    const auto& shard = sub.getShard();
    const size_t nVerts = shard.getVertexCount();
    const size_t globalSize = fes.getSize();

    for (size_t i = 0; i < nVerts; ++i)
    {
      const auto& dofs = fes.getDOFs(0, static_cast<Index>(i));
      for (const Index d : dofs)
      {
        EXPECT_LT(d, globalSize)
          << "P1 DOF " << d << " out of range for local vertex " << i
          << " (orphan-vertex regression, 8x8 triangle mesh).";
      }
    }

    // Also verify owned DOFs are globally unique.
    std::vector<Index> ownedDofs;
    for (size_t i = 0; i < nVerts; ++i)
    {
      if (shard.isOwned(0, i))
      {
        const auto& dofs = fes.getDOFs(0, static_cast<Index>(i));
        for (const Index d : dofs)
          ownedDofs.push_back(d);
      }
    }

    const size_t globalBoundaryVerts = sub.getVertexCount();

    std::vector<std::vector<Index>> allDofs;
    boost::mpi::gather(world, ownedDofs, allDofs, 0);

    if (world.rank() == 0)
    {
      std::vector<Index> combined;
      for (const auto& v : allDofs)
        combined.insert(combined.end(), v.begin(), v.end());

      std::sort(combined.begin(), combined.end());
      const size_t uniqueCount = static_cast<size_t>(
        std::unique(combined.begin(), combined.end()) - combined.begin());

      EXPECT_EQ(uniqueCount, combined.size())
        << "Duplicate owned P1 DOFs on 8x8 boundary SubMesh (orphan-vertex regression).";
      EXPECT_EQ(combined.size(), globalBoundaryVerts)
        << "Owned P1 DOF count does not match global boundary vertex count.";
    }
  }
}

// ---------------------------------------------------------------------------
// main() — initializes MPI environment used by all tests.
// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  g_env = &env;
  g_world = &world;

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
