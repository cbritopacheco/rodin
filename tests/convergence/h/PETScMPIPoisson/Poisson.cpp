/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Distributed PETSc P1/P2 Poisson convergence and global norms. */

#include <cassert>
#include <cmath>
#include <cstdint>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <petsc.h>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/vector.hpp>

#include "Convergence.h"
#include "Rodin/Assembly.h"
#include "Rodin/PETSc.h"
#include "Rodin/MPI/Context/MPI.h"
#include "Rodin/MPI/Geometry/Mesh.h"
#include "Rodin/MPI/Geometry/Sharder.h"
#include "Rodin/MPI/Geometry/SubMesh.h"
#include "Rodin/MPI/Variational/P1.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#include "Rodin/Geometry/BalancedCompactPartitioner.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::PETScMPIPoisson
{
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;

  Mesh<Context::MPI> distribute(
    const Context::MPI& context, Polytope::Type geometry, size_t level)
  {
    Sharder<Context::MPI> sharder(context);
    if (world->rank() == 0)
    {
      auto mesh = UniformGrid(geometry).makeMesh(level);
      const size_t dim = mesh.getDimension();
      mesh.getConnectivity().compute(dim, dim);
      mesh.getConnectivity().compute(dim, 0);
      mesh.getConnectivity().compute(dim, dim - 1);
      mesh.getConnectivity().compute(dim - 1, dim);
      mesh.getConnectivity().compute(dim - 1, 0);
      BalancedCompactPartitioner partitioner(mesh);
      partitioner.partition(static_cast<size_t>(world->size()));
      sharder.shard(partitioner);
      sharder.scatter(0);
    }
    return sharder.gather(0);
  }

  template <class GF, class Exact, class Gradient>
  ErrorNorms globalError(const Mesh<Context::MPI>& mesh, const GF& uh,
    const Exact& exact, const Gradient& exactGradient)
  {
    const auto& shard = mesh.getShard();
    const size_t dim = shard.getDimension();
    const auto discreteGradient = Grad(uh);
    Real localL2 = 0;
    Real localH1 = 0;
    for (Index i = 0; i < shard.getCellCount(); ++i)
    {
      if (!shard.isOwned(dim, i))
        continue;
      const auto cell = shard.getCell(i);
      const auto& qf = QF::PolytopeQuadratureFormula::get(12, cell->getGeometry());
      const auto& quadrature = cell->getQuadrature(qf);
      for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
      {
        const auto& p = quadrature.getPoint(qp);
        const IntegrationPoint ip(p, &qf, qp);
        const Real weight = qf.getWeight(qp) * p.getDistortion();
        const Real valueError = uh(ip) - exact(p);
        localL2 += weight * valueError * valueError;
        const auto gradientError = discreteGradient(ip) - exactGradient(p);
        for (size_t d = 0; d < dim; ++d)
          localH1 += weight * gradientError(d) * gradientError(d);
      }
    }
    const Real totalL2 = boost::mpi::all_reduce(*world, localL2, std::plus<Real>());
    const Real totalH1 = boost::mpi::all_reduce(*world, localH1, std::plus<Real>());
    return {std::sqrt(totalL2), std::sqrt(totalH1)};
  }

  class PETScMPIPoissonTest : public ::testing::TestWithParam<Polytope::Type> {};

  TEST_P(PETScMPIPoissonTest, P2BoundaryConstraintsReachDOFOwners)
  {
    if (world->size() > 4)
      GTEST_SKIP() << "Test designed for at most four ranks.";
    const auto geometry = GetParam();
    const size_t level = geometry == Polytope::Type::Tetrahedron ? 9 : 5;
    Context::MPI context(*environment, *world);
    auto mesh = distribute(context, geometry, level);
    H1<2, Real, Mesh<Context::MPI>> space(
      std::integral_constant<size_t, 2>{}, mesh);
    PETSc::Variational::TrialFunction u(space);
    auto dbc = DirichletBC(u, RealFunction(1));
    dbc.assemble();
    const auto& dofs = std::get<DirichletBCBase<Real>::ValueDOFs>(dbc.getDOFs());
    std::vector<Index> local;
    for (const auto& [index, value] : dofs)
      local.push_back(index);
    std::vector<std::vector<Index>> gathered;
    boost::mpi::all_gather(*world, local, gathered);
    std::set<Index> boundary;
    for (const auto& indices : gathered)
      boundary.insert(indices.begin(), indices.end());
    ASSERT_FALSE(boundary.empty());
    Index begin = 0, end = 0;
    space.getOwnershipRange(begin, end);
    for (const Index index : boundary)
    {
      if (index >= begin && index < end)
      {
        SCOPED_TRACE(::testing::Message() << "rank " << world->rank()
          << ", geometry " << UniformGrid::getGeometryName(geometry)
          << ", boundary DOF " << index);
        EXPECT_TRUE(dofs.contains(index));
      }
    }
  }

  /**
   * On every cell geometry, ranks need P2 rows for owned DOFs and DOFs of
   * owned cells, including off-process slaves, but not unrelated rows.
   */
  TEST_P(PETScMPIPoissonTest, P2IdentificationRowsReachRequiredDOFs)
  {
    if (world->size() > 4)
      GTEST_SKIP() << "Test designed for at most four ranks.";
    const auto geometry = GetParam();
    const size_t level = geometry == Polytope::Type::Tetrahedron ? 9 : 5;
    Context::MPI context(*environment, *world);
    auto mesh = distribute(context, geometry, level);
    H1<2, Real, Mesh<Context::MPI>> space(
      std::integral_constant<size_t, 2>{}, mesh);
    PETSc::Variational::TrialFunction u(space);
    PETSc::Variational::TrialFunction master(space);
    auto dbc = DirichletBC(u, -master);
    dbc.assemble();
    const auto& rows =
      std::get<DirichletBCBase<Real>::IdentifiedDOFs>(dbc.getDOFs());
    std::vector<Index> local;
    for (const auto& [slave, row] : rows)
      local.push_back(slave);
    std::vector<std::vector<Index>> gathered;
    boost::mpi::all_gather(*world, local, gathered);
    std::set<Index> all;
    for (const auto& indices : gathered)
      all.insert(indices.begin(), indices.end());
    ASSERT_FALSE(all.empty());
    std::set<Index> required;
    Index begin, end;
    space.getOwnershipRange(begin, end);
    for (Index i = begin; i < end; ++i)
      required.insert(i);
    for (auto cell = mesh.getCell(); cell; ++cell)
      if (mesh.getShard().isOwned(mesh.getDimension(), cell->getIndex()))
        for (Index dof : space.getDOFs(mesh.getDimension(), cell->getIndex()))
          required.insert(dof);
    for (const Index slave : all)
      EXPECT_EQ(rows.contains(slave), required.contains(slave)) << "rank=" << world->rank()
        << " geometry=" << UniformGrid::getGeometryName(geometry)
        << " slave=" << slave;
  }

  /** Exact P2 polynomial, with either prescribed or affine-identified trace. */
  void checkP2QuadraticPatch(Polytope::Type geometry, bool affine, bool extracted = false)
  {
    if (world->size() > 4)
      GTEST_SKIP() << "Test designed for at most four ranks.";
    const size_t dim = Polytope::Traits(geometry).getDimension();
    const size_t level = geometry == Polytope::Type::Tetrahedron ? 9 : 5;
    const RealFunction exact([dim](const Point& p)
    {
      Real value = 1;
      for (size_t d = 0; d < dim; ++d)
        value += p(d) * p(d);
      return value;
    });
    const RealFunction source(-2 * Real(dim));
    const RealFunction defect([dim](const Point& p)
    {
      Real value = 1;
      for (size_t d = 0; d < dim; ++d)
        value += p(d) * p(d);
      return 2 * value;
    });
    const VectorFunction gradient(dim, [dim](const Point& p)
    {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      for (size_t d = 0; d < dim; ++d)
        value(d) = 2 * p(d);
      return value;
    });
    Context::MPI context(*environment, *world);
    auto parent = distribute(context, geometry, level);
    Optional<SubMesh<Context::MPI>> sub;
    if (extracted)
    {
      SubMesh<Context::MPI>::Builder builder;
      builder.initialize(parent);
      for (auto cell = parent.getCell(); cell; ++cell)
        if (parent.getShard().isOwned(dim, cell->getIndex()))
          builder.include(dim, cell->getIndex());
      sub.emplace(builder.finalize());
    }
    const Mesh<Context::MPI>& mesh = extracted ? static_cast<const Mesh<Context::MPI>&>(*sub) : parent;
    H1<2, Real, Mesh<Context::MPI>> space(
      std::integral_constant<size_t, 2>{}, mesh);
    PETSc::Variational::TrialFunction u(space);
    PETSc::Variational::TestFunction v(space);
    auto stiffness = Integral(Grad(u), Grad(v));
    auto load = Integral(source, v);
    stiffness.setOrder(12);
    load.setOrder(12);
    Problem problem(u, v);
    if (affine)
      problem = stiffness - load + DirichletBC(u, -u, defect);
    else
      problem = stiffness - load + DirichletBC(u, exact);
    PETSc::Solver::CG solver(problem);
    solver.setTolerances(1e-13, 1e-14, 1e5, 20000);
    solver.solve();
    const auto error = globalError(mesh, u.getSolution(), exact, gradient);
    EXPECT_LT(error.getL2(), 1e-8);
    EXPECT_LT(error.getH1Seminorm(), 1e-8);
  }

  TEST_P(PETScMPIPoissonTest, P2QuadraticPatchWithNonzeroDirichletData)
  {
    checkP2QuadraticPatch(GetParam(), false);
  }

  TEST_P(PETScMPIPoissonTest, P2AffineIdentificationQuadraticPatch)
  {
    checkP2QuadraticPatch(GetParam(), true);
  }

  TEST_P(PETScMPIPoissonTest, P2SubMeshQuadraticPatch)
  {
    checkP2QuadraticPatch(GetParam(), false, true);
    checkP2QuadraticPatch(GetParam(), true, true);
  }

  TEST_P(PETScMPIPoissonTest, P2OptimalRates)
  {
    if (world->size() > 4)
      GTEST_SKIP() << "Test designed for at most four ranks.";
    const auto geometry = GetParam();
    const size_t dim = Polytope::Traits(geometry).getDimension();
    const Real pi = Math::Constants::pi();
    const RealFunction exact([dim, pi](const Point& p)
    {
      Real value = 1;
      for (size_t d = 0; d < dim; ++d)
        value *= std::sin(pi * p(d));
      return value;
    });
    const RealFunction source([dim, pi, &exact](const Point& p)
    {
      return Real(dim) * pi * pi * exact(p);
    });
    const VectorFunction gradient(dim, [dim, pi](const Point& p)
    {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      for (size_t d = 0; d < dim; ++d)
      {
        value(d) = pi * std::cos(pi * p(d));
        for (size_t j = 0; j < dim; ++j)
          if (j != d)
            value(d) *= std::sin(pi * p(j));
      }
      return value;
    });
    Context::MPI context(*environment, *world);
    ErrorHistory history;
    const std::vector<size_t> levels =
      geometry == Polytope::Type::Tetrahedron
      ? std::vector<size_t>{7, 9, 11}
      : std::vector<size_t>{5, 7, 9};
    for (const size_t level : levels)
    {
      auto mesh = distribute(context, geometry, level);
      H1<2, Real, Mesh<Context::MPI>> space(
        std::integral_constant<size_t, 2>{}, mesh);
      PETSc::Variational::TrialFunction u(space);
      PETSc::Variational::TestFunction v(space);
      auto stiffness = Integral(Grad(u), Grad(v));
      auto load = Integral(source, v);
      stiffness.setOrder(12);
      load.setOrder(12);
      Problem problem(u, v);
      problem = stiffness - load + DirichletBC(u, exact);
      PETSc::Solver::CG solver(problem);
      solver.setTolerances(1e-12, 1e-14, 1e5, 20000);
      solver.solve();
      EXPECT_TRUE(std::isfinite(solver.getError()));
      history.append(Real(1) / Real(level - 1),
        globalError(mesh, u.getSolution(), exact, gradient));
    }
    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
      const auto rate = history.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message() << "L2 " << coarse.getL2()
        << " -> " << fine.getL2() << ", H1 " << coarse.getH1Seminorm()
        << " -> " << fine.getH1Seminorm() << ", rates " << rate.getL2()
        << ", " << rate.getH1Seminorm());
      EXPECT_GT(rate.getL2(), 2.4);
      EXPECT_LT(rate.getL2(), 3.6);
      EXPECT_GT(rate.getH1Seminorm(), 1.5);
      EXPECT_LT(rate.getH1Seminorm(), 2.5);
    }
  }

  TEST_P(PETScMPIPoissonTest, P1OptimalRates)
  {
    if (world->size() > 4)
      GTEST_SKIP() << "Test designed for at most four ranks.";
    const auto geometry = GetParam();
    const size_t dim = Polytope::Traits(geometry).getDimension();
    const Real pi = Math::Constants::pi();
    const RealFunction exact([dim, pi](const Point& p)
    {
      Real value = 1;
      for (size_t i = 0; i < dim; ++i)
        value *= std::sin(pi * p(i));
      return value;
    });
    const RealFunction source([dim, pi](const Point& p)
    {
      Real value = Real(dim) * pi * pi;
      for (size_t i = 0; i < dim; ++i)
        value *= std::sin(pi * p(i));
      return value;
    });
    const VectorFunction gradient(dim, [dim, pi](const Point& p)
    {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      for (size_t i = 0; i < dim; ++i)
      {
        value(i) = pi * std::cos(pi * p(i));
        for (size_t j = 0; j < dim; ++j)
          if (j != i)
            value(i) *= std::sin(pi * p(j));
      }
      return value;
    });
    Context::MPI context(*environment, *world);
    ErrorHistory history;
    for (size_t level : {size_t(5), size_t(9), size_t(17)})
    {
      auto mesh = distribute(context, geometry, level);
      P1<Real, Mesh<Context::MPI>> space(mesh);
      PETSc::Variational::TrialFunction u(space);
      PETSc::Variational::TestFunction v(space);
      auto stiffness = Integral(Grad(u), Grad(v));
      auto load = Integral(source, v);
      stiffness.setOrder(12);
      load.setOrder(12);
      Problem problem(u, v);
      problem = stiffness - load + DirichletBC(u, exact);
      PETSc::Solver::CG solver(problem);
      solver.setTolerances(1e-12, 1e-14, 1e5, 20000);
      solver.solve();
      EXPECT_TRUE(std::isfinite(solver.getError()));
      history.append(Real(1) / Real(level - 1),
        globalError(mesh, u.getSolution(), exact, gradient));
    }
    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
      const auto rate = history.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message() << "L2 " << coarse.getL2()
        << " -> " << fine.getL2() << ", H1 " << coarse.getH1Seminorm()
        << " -> " << fine.getH1Seminorm() << ", rates " << rate.getL2()
        << ", " << rate.getH1Seminorm());
      EXPECT_GT(rate.getL2(), 1.5);
      EXPECT_LT(rate.getL2(), 2.5);
      EXPECT_GT(rate.getH1Seminorm(), 0.7);
      EXPECT_LT(rate.getH1Seminorm(), 1.5);
    }
  }

  /**
   * Manufactured identification rates for the two previously affected
   * boundary paths.  The affine case shifts the exact solution by one:
   * the trace condition u = -u + 2 then prescribes u = 1 on the boundary.
   */
  template <size_t Order, bool Affine>
  void checkIdentificationRates(Polytope::Type geometry)
  {
    if (world->size() > 4)
      GTEST_SKIP() << "Test designed for at most four ranks.";
    const size_t dim = Polytope::Traits(geometry).getDimension();
    const Real pi = Math::Constants::pi();
    const Real offset = Affine ? Real(1) : Real(0);
    const RealFunction exact([dim, pi, offset](const Point& p)
    {
      Real value = 1;
      for (size_t d = 0; d < dim; ++d)
        value *= std::sin(pi * p(d));
      return offset + value;
    });
    const RealFunction source([dim, pi](const Point& p)
    {
      Real value = Real(dim) * pi * pi;
      for (size_t d = 0; d < dim; ++d)
        value *= std::sin(pi * p(d));
      return value;
    });
    const VectorFunction gradient(dim, [dim, pi](const Point& p)
    {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      for (size_t d = 0; d < dim; ++d)
      {
        value(d) = pi * std::cos(pi * p(d));
        for (size_t j = 0; j < dim; ++j)
          if (j != d)
            value(d) *= std::sin(pi * p(j));
      }
      return value;
    });
    const std::vector<size_t> levels = Order == 1
      ? std::vector<size_t>{5, 9, 17}
      : geometry == Polytope::Type::Tetrahedron
        ? std::vector<size_t>{7, 9, 11}
        : std::vector<size_t>{5, 7, 9};
    Context::MPI context(*environment, *world);
    ErrorHistory history;
    for (const size_t level : levels)
    {
      auto mesh = distribute(context, geometry, level);
      auto solveLevel = [&](auto& space)
      {
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
        auto stiffness = Integral(Grad(u), Grad(v));
        auto load = Integral(source, v);
        stiffness.setOrder(12);
        load.setOrder(12);
        Problem problem(u, v);
        if constexpr (Affine)
          problem = stiffness - load
            + DirichletBC(u, -u, RealFunction(2));
        else
          problem = stiffness - load + DirichletBC(u, -u);
        PETSc::Solver::CG solver(problem);
        solver.setTolerances(1e-12, 1e-14, 1e5, 20000);
        solver.solve();
        EXPECT_TRUE(std::isfinite(solver.getError()));
        history.append(Real(1) / Real(level - 1),
          globalError(mesh, u.getSolution(), exact, gradient));
      };
      if constexpr (Order == 1)
      {
        P1<Real, Mesh<Context::MPI>> space(mesh);
        solveLevel(space);
      }
      else
      {
        H1<2, Real, Mesh<Context::MPI>> space(
          std::integral_constant<size_t, 2>{}, mesh);
        solveLevel(space);
      }
    }
    ASSERT_EQ(history.getSize(), 3u);
    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
      const auto rate = history.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message()
        << "geometry=" << UniformGrid::getGeometryName(geometry)
        << " order=" << Order << " affine=" << Affine
        << " L2 " << coarse.getL2() << " -> " << fine.getL2()
        << " H1 " << coarse.getH1Seminorm() << " -> "
        << fine.getH1Seminorm() << " rates " << rate.getL2()
        << ", " << rate.getH1Seminorm());
      if constexpr (Order == 1)
      {
        EXPECT_GT(rate.getL2(), 1.5);
        EXPECT_LT(rate.getL2(), 2.5);
        EXPECT_GT(rate.getH1Seminorm(), 0.7);
        EXPECT_LT(rate.getH1Seminorm(), 1.5);
      }
      else
      {
        EXPECT_GT(rate.getL2(), 2.4);
        EXPECT_LT(rate.getL2(), 3.6);
        EXPECT_GT(rate.getH1Seminorm(), 1.5);
        EXPECT_LT(rate.getH1Seminorm(), 2.5);
      }
    }
  }

  TEST_P(PETScMPIPoissonTest, P1HomogeneousIdentificationOptimalRates)
  {
    checkIdentificationRates<1, false>(GetParam());
  }

  TEST_P(PETScMPIPoissonTest, P2AffineIdentificationOptimalRates)
  {
    checkIdentificationRates<2, true>(GetParam());
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScMPIPoissonTest,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const ::testing::TestParamInfo<Polytope::Type>& info)
    {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}

int main(int argc, char** argv)
{
  [[maybe_unused]] const PetscErrorCode ierr =
    PetscInitialize(&argc, &argv, nullptr, nullptr);
  assert(ierr == PETSC_SUCCESS);
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator comm;
  Rodin::Tests::Convergence::H::PETScMPIPoisson::environment = &env;
  Rodin::Tests::Convergence::H::PETScMPIPoisson::world = &comm;
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  [[maybe_unused]] const PetscErrorCode finalizeError = PetscFinalize();
  assert(finalizeError == PETSC_SUCCESS);
  return result;
}
