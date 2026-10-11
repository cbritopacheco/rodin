/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * PETSc-backed H1 restriction from distributed parent meshes to full,
 * boundary, sparse and nested SubMeshes. Orders 1-6 and all seven parent
 * geometries use degree-matched reference polynomial oracles on affine and
 * exact quadratic domains. Real and
 * complex PETSc configurations each exercise scalar, vector and matrix fields.
 * Every held coefficient and physical sample is checked independently;
 * correspondence is determined exclusively by logical ancestry.
 * Mesh/space construction, interpolation and ghost-state updates are
 * collective. Synchronized field reads and metadata queries are also
 * exercised on rank zero alone; barriers belong only to the test protocol.
 */
#include <cmath>
#include <cassert>
#include <type_traits>
#include <utility>
#include <vector>
#include <gtest/gtest.h>
#include <petsc.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/BalancedCompactPartitioner.h>
#include <Rodin/MPI/Context/MPI.h>
#include <Rodin/MPI/Geometry/Sharder.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Geometry/SubMesh.h>
#include <Rodin/MPI/Variational/H1.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>
#include "../../../convergence/CurvedGeometry.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
static boost::mpi::environment* g_env = nullptr;
static boost::mpi::communicator* g_world = nullptr;
#ifdef PETSC_USE_COMPLEX
using BackendScalar = Complex;
#else
using BackendScalar = Real;
#endif
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
    Polytope::Type type, std::initializer_list<size_t> shape, bool curved)
  {
    const auto& comm = ctx.getCommunicator();
    Sharder<Context::MPI> sharder(ctx);
    if (comm.rank() == 0)
    {
      auto localMesh = makeShardableMesh(type, shape);
      if (curved)
      {
        // Install before distribution: curved maps must survive both shard
        // transport and subsequent full/sparse/boundary/nested extraction.
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

}

namespace Rodin::Tests::Unit
{
  TEST(PETScMPISubMeshH1, ExplicitRestriction_AllOrdersBackendRangesAllGeometries)
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
            GridFunction<Space, ::Vec> source(sourceSpace), restricted(targetSpace),
              oracle(targetSpace), rejected(targetSpace);
            const auto scalar = [curved](const Point& point, size_t component) {
              Math::SpatialPoint x = point.getPhysicalCoordinates();
              if (curved)
              {
                // Independent analytic inverse of the prescribed global map.
                // This defines field data, not entity/DOF correspondence.
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
            rejected.project(data);
            rejected *= Scalar(0);
            // Interpolation and ghost updates change distributed state collectively.
            restricted.project(source);
            const auto checkValue = [&](const Point& point) {
              const auto actual = std::as_const(restricted)(point),
                         expected = data(point), wrong = std::as_const(rejected)(point);
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
            if (world.rank() == 0)
            {
              // Peers do not read fields here: synchronized held-DOF access,
              // cached metadata and point evaluation must remain local operations.
              EXPECT_EQ(&restricted.getFiniteElementSpace(), &targetSpace);
              EXPECT_EQ(restricted.getSize(), targetSpace.getSize());
              for (auto cell = mesh.getCell(); cell; ++cell)
              {
                const IndexArray dofs =
                  targetSpace.getDOFs(mesh.getDimension(), cell->getIndex());
                for (Index local = 0; local < static_cast<size_t>(dofs.size()); ++local)
                {
                  EXPECT_LT(std::abs(std::as_const(restricted)[dofs(local)] -
                              std::as_const(oracle)[dofs(local)]),
                    PolynomialTolerance);
                }
                const Polytope::Traits traits(cell->getGeometry());
                checkValue(Point(*cell, traits.getCentroid()));
              }
            }
            // Test-protocol synchronization, outside the local read operation.
            world.barrier();
            EXPECT_EQ(&restricted.getFiniteElementSpace(), &targetSpace);
            EXPECT_EQ(restricted.getSize(), targetSpace.getSize());
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
                EXPECT_LT(std::abs(std::as_const(restricted)[dofs(local)] -
                            std::as_const(oracle)[dofs(local)]),
                  PolynomialTolerance);
              }
              const Polytope::Traits traits(cell->getGeometry());
              const auto sample = [&](const Math::SpatialPoint& coordinates) {
                checkValue(Point(*cell, coordinates));
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
            range.template operator()<K, BackendScalar>();
            range.template operator()<K, Math::SpatialVector<BackendScalar>>();
            range.template operator()<K, Math::SpatialMatrix<BackendScalar>>();
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

}

int main(int argc, char** argv)
{
  boost::mpi::environment environment(argc, argv);
  boost::mpi::communicator world;
  g_env = &environment;
  g_world = &world;
  PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
  if (ierr != PETSC_SUCCESS)
    return static_cast<int>(ierr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  ierr = PetscFinalize();
  return result ? result : static_cast<int>(ierr);
}
