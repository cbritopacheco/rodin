/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Matrix finite element, tensor operator, and backend regressions. */
#include <gtest/gtest.h>
#include <boost/mpi.hpp>
#include <cstdio>
#include <sstream>
#include "Rodin/Geometry.h"
#include "Rodin/Geometry/BalancedCompactPartitioner.h"
#include "Rodin/MPI/Context/MPI.h"
#include "Rodin/MPI/Geometry/Sharder.h"
#include "Rodin/MPI/Variational/P0.h"
#include "Rodin/MPI/Variational/P0g.h"
#include "Rodin/MPI/Variational/P1.h"
#include "Rodin/MPI/Variational/H1.h"
#include "Rodin/Variational.h"
#ifdef RODIN_TEST_WITH_PETSC
#include "Rodin/PETSc.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
#ifdef RODIN_TEST_WITH_PETSC
using BackendScalar = PetscScalar;
#else
using BackendScalar = Real;
#endif
using BackendMatrix = Math::SpatialMatrix<BackendScalar>;
static boost::mpi::environment* environment;
static boost::mpi::communicator* world;

namespace
{
  std::vector<std::pair<size_t, size_t>> matrixShapes()
  {
    std::vector<std::pair<size_t, size_t>> shapes;
    for (size_t rows = 1; rows <= 3; ++rows)
      for (size_t cols = 1; cols <= 3; ++cols)
        shapes.emplace_back(rows, cols);
    return shapes;
  }

#ifdef RODIN_TEST_WITH_PETSC
  template <class FES>
  void checkPetscSolve(const FES& fes, bool diffusion)
  {
    auto exact =
      MatrixFunction(fes.getRows(), fes.getColumns(), [&](const Geometry::Point& p) {
        BackendMatrix value(fes.getRows(), fes.getColumns());
        for (size_t r = 0; r < fes.getRows(); ++r)
          for (size_t c = 0; c < fes.getColumns(); ++c)
            value(r, c) = 1 + 7 * r + c + (diffusion ? p.x() : Real(0));
        return value;
      });
    exact.setOrder(diffusion ? 1 : 0);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);
    auto mass = Integral(u, v);
    auto load = Integral(exact, v);
    auto stiffness = Integral(Grad(u), Grad(v));
    // Collapsed pyramid bases require more than the polynomial degree rule.
    mass.setOrder(12);
    load.setOrder(12);
    stiffness.setOrder(12);
    Problem problem(u, v);
    if (diffusion)
      problem = stiffness + mass - load + DirichletBC(u, exact);
    else
      problem = mass - load;
    PETSc::Solver::CG solver(problem);
    solver.setTolerances(1e-11, 1e-13, 1e5, 1000);
    solver.solve();
    for (auto it = fes.getMesh().getCell(); it; ++it)
    {
      Geometry::Point point(*it, Polytope::Traits(it->getGeometry()).getCentroid());
      EXPECT_LE((u.getSolution().getValue(point) - exact.getValue(point)).norm(), 1e-7);
    }
  }
#endif
  template <class MatrixSpace, class ScalarSpace>
  void checkDistributed(const MatrixSpace& matrix, const ScalarSpace& scalar)
  {
    const Index components = matrix.getVectorDimension();
    ASSERT_EQ(matrix.getSize(), scalar.getSize() * components);
    ASSERT_EQ(matrix.getShard().getSize(), scalar.getShard().getSize() * components);
    Index begin, end, scalarBegin, scalarEnd;
    matrix.getOwnershipRange(begin, end);
    scalar.getOwnershipRange(scalarBegin, scalarEnd);
    EXPECT_EQ(begin, scalarBegin * components);
    EXPECT_EQ(end, scalarEnd * components);
    for (size_t local = 0; local < matrix.getShard().getSize(); ++local)
    {
      const Index global = matrix.getGlobalIndex(local);
      EXPECT_EQ(global,
        scalar.getGlobalIndex(local / components) * components + local % components);
      ASSERT_TRUE(matrix.getLocalIndex(global));
      EXPECT_EQ(*matrix.getLocalIndex(global), local);
    }
    EXPECT_FALSE(matrix.getLocalIndex(matrix.getSize()));
#ifdef RODIN_TEST_WITH_PETSC
    BackendMatrix value(matrix.getRows(), matrix.getColumns());
    for (size_t r = 0; r < matrix.getRows(); ++r)
      for (size_t c = 0; c < matrix.getColumns(); ++c)
        value(r, c) = 1 + 10 * r + c;
    GridFunction<MatrixSpace, ::Vec> field(matrix);
#if defined(PETSC_USE_COMPLEX)
    for (size_t r = 0; r < matrix.getRows(); ++r)
      for (size_t c = 0; c < matrix.getColumns(); ++c)
        value(r, c) += PetscCMPLX(0, 1 + r + c);
#endif
    field.project([&](const Geometry::Point&) { return value; });
    PETSc::Variational::TrialFunction petscTrial(matrix);
    PETSc::Variational::TestFunction petscTest(matrix);
    BilinearForm globalMass(petscTrial, petscTest);
    auto massIntegrator = Integral(petscTrial, petscTest);
    massIntegrator.setOrder(12);
    globalMass = massIntegrator;
    globalMass.assemble();
    PetscInt globalRows, globalCols;
    MatGetSize(globalMass.getOperator(), &globalRows, &globalCols);
    EXPECT_EQ(globalRows, matrix.getSize());
    EXPECT_EQ(globalCols, matrix.getSize());
    LinearForm load(petscTest);
    auto loadIntegrator = Integral(MatrixFunction(value), petscTest);
    loadIntegrator.setOrder(12);
    load = loadIntegrator;
    load.assemble();
    ::Vec product = nullptr;
    ASSERT_EQ(VecDuplicate(load.getVector(), &product), PETSC_SUCCESS);
    ASSERT_EQ(MatMult(globalMass.getOperator(), field.getData(), product), PETSC_SUCCESS);
    ASSERT_EQ(VecAXPY(product, -1, load.getVector()), PETSC_SUCCESS);
    PetscReal residual = 0;
    ASSERT_EQ(VecNorm(product, NORM_INFINITY, &residual), PETSC_SUCCESS);
    EXPECT_NEAR(residual, 0, 1e-9);
    ASSERT_EQ(VecDestroy(&product), PETSC_SUCCESS);
    const std::string filename =
      boost::filesystem::unique_path("/tmp/rodin_matrix_petsc_%%%%-%%%%-%%%%.h5")
        .string();
    IO::GridFunctionPrinter<IO::FileFormat::HDF5, MatrixSpace, ::Vec>(field).print(
      filename);
    GridFunction<MatrixSpace, ::Vec> loaded(matrix);
    IO::GridFunctionLoader<IO::FileFormat::HDF5, MatrixSpace, ::Vec>(loaded).load(
      filename);
    PetscReal loadedNorm = 0;
    ASSERT_EQ(VecAXPY(loaded.getData(), -1, field.getData()), PETSC_SUCCESS);
    ASSERT_EQ(VecNorm(loaded.getData(), NORM_INFINITY, &loadedNorm), PETSC_SUCCESS);
    EXPECT_NEAR(loadedNorm, 0, 1e-12);
    std::remove(filename.c_str());
    if constexpr (!std::is_same_v<typename MatrixSpace::ElementType,
                    P0gElement<typename MatrixSpace::RangeType>>)
    {
      std::stringstream stream;
      IO::GridFunctionPrinter<IO::FileFormat::MFEM, MatrixSpace, ::Vec>(field).print(
        stream);
      EXPECT_NE(
        stream.str().find("VDim: " + std::to_string(components)), std::string::npos);
    }

#endif
    TrialFunction trial(matrix);
    TestFunction test(matrix);
    auto mass = Integral(trial, test);
    for (auto it = matrix.getMesh().getCell(); it; ++it)
    {
      const auto& dofs = matrix.getDOFs(it->getDimension(), it->getIndex());
      const auto& localDOFs =
        matrix.getShard().getDOFs(it->getDimension(), it->getIndex());
      for (size_t a = 0; a < static_cast<size_t>(dofs.size()); ++a)
      {
        EXPECT_EQ(dofs[a], matrix.getGlobalIndex(localDOFs[a]));
        EXPECT_EQ(
          dofs[a], matrix.getGlobalIndex({it->getDimension(), it->getIndex()}, a));
      }
      mass.setPolytope(*it);
      EXPECT_GT(std::abs(mass.integrate(0, 0)), 0);
#ifdef RODIN_TEST_WITH_PETSC
      Geometry::Point point(*it, Polytope::Traits(it->getGeometry()).getCentroid());
      const auto actual = field.getValue(point);
      for (size_t r = 0; r < matrix.getRows(); ++r)
        for (size_t c = 0; c < matrix.getColumns(); ++c)
          EXPECT_NEAR(std::abs(actual(r, c) - value(r, c)), 0, 1e-10);
#endif
    }
#ifdef RODIN_TEST_WITH_PETSC
    checkPetscSolve(matrix,
      !std::is_same_v<typename MatrixSpace::ElementType,
        P0Element<typename MatrixSpace::RangeType>> &&
        !std::is_same_v<typename MatrixSpace::ElementType,
          P0gElement<typename MatrixSpace::RangeType>>);
#endif
  }
}

TEST(DistributedMatrixRange, AllSpaces)
{
  Context::MPI context(*environment, *world);
  for (auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Pyramid, Polytope::Type::Wedge, Polytope::Type::Hexahedron})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    Sharder<Context::MPI> sharder(context);
    if (world->rank() == 0)
    {
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      auto mesh = dimension == 1 ? LocalMesh::UniformGrid(geometry, {9})
        : dimension == 2         ? LocalMesh::UniformGrid(geometry, {4, 4})
                                 : LocalMesh::UniformGrid(geometry, {3, 3, 3});
      const size_t D = mesh.getDimension();
      for (size_t d = 1; d <= D; ++d)
        for (size_t lower = 0; lower < d; ++lower)
        {
          mesh.getConnectivity().compute(d, lower);
          mesh.getConnectivity().compute(lower, d);
        }
      mesh.getConnectivity().compute(D, D);
      mesh.getConnectivity().compute(D, 0);
      mesh.getConnectivity().compute(D, D - 1);
      mesh.getConnectivity().compute(D - 1, D);
      mesh.getConnectivity().compute(D - 1, 0);
      BalancedCompactPartitioner partitioner(mesh);
      partitioner.partition(world->size());
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    auto mesh = sharder.gather(0);
    auto& shard = mesh.getShard();
    for (size_t d = 1; d <= shard.getDimension(); ++d)
      for (size_t lower = 0; lower < d; ++lower)
        shard.getConnectivity().compute(d, lower);
    for (const auto shape : matrixShapes())
    {
      const auto [rows, cols] = shape;
      P0<BackendScalar, Geometry::Mesh<Context::MPI>> p0(mesh);
      P0<BackendMatrix, Geometry::Mesh<Context::MPI>> mp0(mesh, rows, cols);
      checkDistributed(mp0, p0);
      P0g<BackendScalar, Geometry::Mesh<Context::MPI>> p0g(mesh);
      P0g<BackendMatrix, Geometry::Mesh<Context::MPI>> mp0g(mesh, rows, cols);
      checkDistributed(mp0g, p0g);
      P1<BackendScalar, Geometry::Mesh<Context::MPI>> p1(mesh);
      P1<BackendMatrix, Geometry::Mesh<Context::MPI>> mp1(mesh, rows, cols);
      checkDistributed(mp1, p1);
      H1<1, BackendScalar, Geometry::Mesh<Context::MPI>> h1(
        std::integral_constant<size_t, 1>{}, mesh);
      H1<1, BackendMatrix, Geometry::Mesh<Context::MPI>> mh1(
        std::integral_constant<size_t, 1>{}, mesh, rows, cols);
      checkDistributed(mh1, h1);
      H1<2, BackendScalar, Geometry::Mesh<Context::MPI>> h2(
        std::integral_constant<size_t, 2>{}, mesh);
      H1<2, BackendMatrix, Geometry::Mesh<Context::MPI>> mh2(
        std::integral_constant<size_t, 2>{}, mesh, rows, cols);
      checkDistributed(mh2, h2);
      H1<3, BackendScalar, Geometry::Mesh<Context::MPI>> h3(
        std::integral_constant<size_t, 3>{}, mesh);
      H1<3, BackendMatrix, Geometry::Mesh<Context::MPI>> mh3(
        std::integral_constant<size_t, 3>{}, mesh, rows, cols);
      checkDistributed(mh3, h3);
    }
  }
}

#ifdef RODIN_TEST_WITH_PETSC
TEST(LocalPetscMatrixRange, AllSpaces)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  for (size_t d = 1; d <= 2; ++d)
    for (size_t lower = 0; lower < d; ++lower)
    {
      mesh.getConnectivity().compute(d, lower);
      mesh.getConnectivity().compute(lower, d);
    }
  checkPetscSolve(P0<BackendMatrix, LocalMesh>(mesh, 2, 3), false);
  checkPetscSolve(P0g<BackendMatrix, LocalMesh>(mesh, 2, 3), false);
  checkPetscSolve(P1<BackendMatrix, LocalMesh>(mesh, 2, 3), true);
  checkPetscSolve(
    H1<1, BackendMatrix>(std::integral_constant<size_t, 1>{}, mesh, 2, 3), true);
  checkPetscSolve(
    H1<2, BackendMatrix>(std::integral_constant<size_t, 2>{}, mesh, 2, 3), true);
  checkPetscSolve(
    H1<3, BackendMatrix>(std::integral_constant<size_t, 3>{}, mesh, 2, 3), true);
}
#endif

int main(int argc, char** argv)
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator comm;
  environment = &env;
  world = &comm;
#ifdef RODIN_TEST_WITH_PETSC
  PetscInitialize(&argc, &argv, nullptr, nullptr);
#endif
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
#ifdef RODIN_TEST_WITH_PETSC
  PetscFinalize();
#endif
  return result;
}
