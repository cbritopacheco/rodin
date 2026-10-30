/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "Rodin/Assembly/Sequential.h"
#ifdef RODIN_USE_OPENMP
#include "Rodin/Assembly/OpenMP.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /// @brief Compares constrained mass solutions with independent component algebra.
  template <class Operator, template <class, class> class Backend>
  void checkConstraintChains()
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
    mesh.getConnectivity().compute(0, 1);
    P1<Math::SpatialVector<Real>> space(mesh, 3);
    TrialFunction u(space);
    TestFunction v(space);
    using System = Math::LinearSystem<Operator, Math::Vector<Real>>;
    using ProblemType = Problem<System, decltype(u), decltype(v)>;
    using Body = typename ProblemType::ProblemBodyType;
    Backend<System, ProblemType> assembler;
    if constexpr (requires { assembler.setThreadCount(size_t(2)); })
      assembler.setThreadCount(2);
    VectorFunction forcing{3., 0., 0.};
    for (const bool zero : {false, true})
    {
      for (const bool affine : {false, true})
      {
        SCOPED_TRACE(zero);
        SCOPED_TRACE(affine);
        Math::SpatialMatrix<Real> transform(3, 3);
        transform.setZero();
        if (!zero)
        {
          transform(0, 1) = 1;
          transform(1, 2) = 1;
          transform(2, 2) = 1;
        }
        const Math::SpatialVector<Real> defect = affine
          ? (zero ? Math::SpatialVector<Real>{7, 8, 9}
                  : Math::SpatialVector<Real>{1, 2, 0})
          : Math::SpatialVector<Real>{0, 0, 0};
        const auto offset =
          VectorFunction(size_t(3), [defect](const Point&) { return defect; });
        Body body = Integral(u, v) - Integral(forcing, v) +
          DirichletBC(u, MatrixFunction(transform) * u, offset);
        Assembly::ProblemAssemblyInput<Body, decltype(u), decltype(v)> input(body, u, v);
        System system;
        assembler.execute(system, input);
        const Math::Matrix<Real> matrix = system.getOperator();
        const Math::Vector<Real> solution = matrix.fullPivLu().solve(system.getVector());
        const Math::SpatialVector<Real> expected = zero
          ? defect
          : (affine ? Math::SpatialVector<Real>{Real(7) / 3, Real(4) / 3, -Real(2) / 3}
                    : Math::SpatialVector<Real>{1, 1, 1});
        const auto cell = mesh.getCell();
        for (size_t local = 0; local < u.getDOFs(*cell); ++local)
        {
          EXPECT_NEAR(solution[space.getGlobalIndex({1, cell->getIndex()}, local)],
            expected[local % 3], 1e-12);
        }
      }
    }
    Math::SpatialMatrix<Real> cycle(3, 3);
    cycle.setZero();
    cycle(0, 1) = cycle(1, 0) = cycle(2, 2) = 1;
    Body cyclic =
      Integral(u, v) - Integral(forcing, v) + DirichletBC(u, MatrixFunction(cycle) * u);
    Assembly::ProblemAssemblyInput<Body, decltype(u), decltype(v)> input(cyclic, u, v);
    System system;
    EXPECT_THROW(assembler.execute(system, input), Alert::Exception);
  }
}

/// @brief Dense sequential assembly resolves affine chains and zero identification rows.
TEST(ConstraintChain, SequentialDense)
{
  checkConstraintChains<Math::Matrix<Real>, Assembly::Sequential>();
}

/// @brief Sparse sequential assembly resolves affine chains and zero identification rows.
TEST(ConstraintChain, SequentialSparse)
{
  checkConstraintChains<Math::SparseMatrix<Real>, Assembly::Sequential>();
}

#ifdef RODIN_USE_OPENMP
/// @brief Dense OpenMP assembly resolves affine chains and zero identification rows.
TEST(ConstraintChain, OpenMPDense)
{
  checkConstraintChains<Math::Matrix<Real>, Assembly::OpenMP>();
}

/// @brief Sparse OpenMP assembly resolves affine chains and zero identification rows.
TEST(ConstraintChain, OpenMPSparse)
{
  checkConstraintChains<Math::SparseMatrix<Real>, Assembly::OpenMP>();
}
#endif

/// @brief Empty named forms do not form a NaN reserve estimate.
TEST(ConstraintChain, EmptyNamedForm)
{
  const auto mesh = LocalMesh::Builder().initialize(1).nodes(0).finalize();
  P1<Real> space(mesh);
  TrialFunction trial(space);
  TestFunction test(space);
  const MassForm mass(trial, test);
  EXPECT_EQ(mass.getOperator().rows(), 0);
  EXPECT_EQ(mass.getOperator().cols(), 0);
}
