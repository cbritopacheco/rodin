/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Distributed PETSc P1 Poisson convergence and global norms. */

#include <cassert>
#include <cmath>
#include <cstdint>
#include <string>

#include <gtest/gtest.h>
#include <petsc.h>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include "Convergence.h"
#include "Rodin/Assembly.h"
#include "Rodin/PETSc.h"
#include "Rodin/MPI/Context/MPI.h"
#include "Rodin/MPI/Geometry/Mesh.h"
#include "Rodin/MPI/Geometry/Sharder.h"
#include "Rodin/MPI/Variational/P1.h"
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
