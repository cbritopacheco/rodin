/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Independent mass-system oracles for PETSc problem ownership and target updates.
 *
 * Partial assembly retains full Dirichlet elimination and raw-form lifting.
 * The unselected stored matrix or vector remains unchanged, including after
 * external scaling.
 */
#include <gtest/gtest.h>
#include <array>
#include <memory>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>
#include "Rodin/Geometry.h"
#include "Rodin/PETSc.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  using Scalar = PetscScalar;

  /// @brief Reads a local matrix without depending on its sparse storage layout.
  Math::Matrix<Scalar> matrix(const PETSc::Math::LinearSystem& system)
  {
    PetscInt rows = 0, cols = 0;
    EXPECT_EQ(MatGetSize(system.getOperator(), &rows, &cols), PETSC_SUCCESS);
    Math::Matrix<Scalar> values(rows, cols);
    for (PetscInt i = 0; i < rows; ++i)
    {
      for (PetscInt j = 0; j < cols; ++j)
      {
        EXPECT_EQ(
          MatGetValues(system.getOperator(), 1, &i, 1, &j, &values(i, j)), PETSC_SUCCESS);
      }
    }
    return values;
  }

  /// @brief Reads an owned sequential vector into an independent value snapshot.
  Math::Vector<Scalar> vector(::Vec source)
  {
    PetscInt size = 0;
    EXPECT_EQ(VecGetSize(source, &size), PETSC_SUCCESS);
    std::vector<PetscInt> indices(static_cast<size_t>(size));
    std::iota(indices.begin(), indices.end(), 0);
    Math::Vector<Scalar> values(size);
    EXPECT_EQ(VecGetValues(source, size, indices.data(), values.data()), PETSC_SUCCESS);
    return values;
  }

  /// @brief Problem copies retain shared systems and refresh a stale body before solving.
  template <class ProblemType, class CheckFields>
  void checkOwnership(ProblemType& original, const Math::Matrix<Scalar>& raw,
    const Math::Vector<Scalar>& load, const Math::Vector<Scalar>& exact,
    const std::vector<Index>& fixed, const CheckFields& checkFields)
  {
    const Eigen::Index size = raw.rows();
    Math::Matrix<Scalar> projection = Math::Matrix<Scalar>::Identity(size, size);
    Math::Vector<Scalar> lift = Math::Vector<Scalar>::Zero(size);
    for (const auto dof : fixed)
    {
      projection(dof, dof) = 0;
      lift[dof] = exact[dof];
    }
    const Math::Matrix<Scalar> identity = Math::Matrix<Scalar>::Identity(size, size);
    const Math::Matrix<Scalar> full =
      projection * raw * projection + identity - projection;
    const Math::Vector<Scalar> fullLoad = projection * (load - raw * lift) + lift;
    EXPECT_LT((raw * exact - load).norm(), 1e-13);
    auto checkFull = [&](const auto& problem) {
      const auto& system = problem.getLinearSystem();
      EXPECT_LT((matrix(system) - full).norm(), 1e-12);
      EXPECT_LT((vector(system.getVector()) - fullLoad).norm(), 1e-12);
    };
    checkFull(std::as_const(original));
    auto body = original.getBody();
    auto verify = [&](auto& problem, const char* stage) {
      SCOPED_TRACE(stage);
      checkFull(std::as_const(problem));
      auto& system = problem.getLinearSystem();
      EXPECT_EQ(system.getOperator(), original.getLinearSystem().getOperator());
      EXPECT_EQ(system.getVector(), original.getLinearSystem().getVector());
      EXPECT_EQ(system.getSolution(), original.getLinearSystem().getSolution());
      ASSERT_EQ(VecSet(system.getVector(), Scalar(7)), PETSC_SUCCESS);
      problem.assemble(AssemblyTarget::LHS);
      EXPECT_LT((matrix(system) - full).norm(), 1e-12);
      EXPECT_LT(
        (vector(system.getVector()) - Math::Vector<Scalar>::Constant(size, 7)).norm(),
        1e-12);
      EXPECT_LT((matrix(std::as_const(original).getLinearSystem()) - full).norm(), 1e-12);
      EXPECT_LT((vector(std::as_const(original).getLinearSystem().getVector()) -
                  Math::Vector<Scalar>::Constant(size, 7))
                  .norm(),
        1e-12);
      problem.assemble();
      ASSERT_EQ(MatScale(system.getOperator(), Scalar(3)), PETSC_SUCCESS);
      problem.assemble(AssemblyTarget::RHS);
      EXPECT_LT((matrix(system) - Scalar(3) * full).norm(), 1e-12);
      EXPECT_LT((vector(system.getVector()) - fullLoad).norm(), 1e-12);
      EXPECT_LT(
        (matrix(std::as_const(original).getLinearSystem()) - Scalar(3) * full).norm(),
        1e-12);
      EXPECT_LT(
        (vector(std::as_const(original).getLinearSystem().getVector()) - fullLoad).norm(),
        1e-12);
      std::unique_ptr<typename ProblemType::ProblemBodyType> refreshed(body.copy());
      problem = *refreshed;
      refreshed.reset();
      PETSc::Solver::GMRES solver(problem);
      solver.setTolerances(1e-12, 1e-14, 1e6, 100);
      ::PC preconditioner = nullptr;
      ASSERT_EQ(KSPGetPC(solver.getHandle(), &preconditioner), PETSC_SUCCESS);
      ASSERT_EQ(PCSetType(preconditioner, PCJACOBI), PETSC_SUCCESS);
      solver.solve();
      EXPECT_LT((vector(system.getSolution()) - exact).norm(), 1e-10);
      checkFields();
      checkFull(std::as_const(problem));
      checkFull(std::as_const(original));
    };
    auto copied = original;
    verify(copied, "copy construction");
    auto moved = std::move(copied);
    verify(moved, "move construction");
    std::unique_ptr<ProblemType> cloned(original.copy());
    verify(*cloned, "polymorphic clone");
    auto assigned = original;
    assigned = original;
    verify(assigned, "copy assignment");
    auto moveSource = original;
    auto moveAssigned = original;
    moveAssigned = std::move(moveSource);
    verify(moveAssigned, "move assignment");
    std::unique_ptr<ProblemType> survivor;
    {
      const ProblemType temporary(original);
      survivor.reset(temporary.copy());
    }
    verify(*survivor, "source destroyed");
  }

  /// @brief Single and mixed scalar fields obey independent closed-form interval mass matrices.
  template <bool Mixed, class FieldScalar = Scalar>
  void checkProblemLifecycle()
  {
    const FieldScalar phase = [] {
      if constexpr (std::is_same_v<FieldScalar, Complex>)
        return FieldScalar(1, -0.5);
      else
        return FieldScalar(1);
    }();
    const auto known = [](FieldScalar value) {
      if constexpr (std::is_same_v<FieldScalar, Complex>)
        return ComplexFunction(value);
      else
        return RealFunction(value);
    };
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    mesh.scale(0.5);
    mesh.getConnectivity().compute(0, 1);
    P1<FieldScalar> space(mesh);
    ASSERT_EQ(space.getSize(), 3u);
    std::array<Index, 3> dofs{};
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
    {
      const Real x = mesh.getVertexCoordinates(vertex->getIndex()).x();
      const size_t node = static_cast<size_t>(2 * x);
      ASSERT_LT(node, dofs.size());
      dofs[node] = space.getDOFs(0, vertex->getIndex())[0];
      if (x == 0)
        mesh.setAttribute(vertex.key(), 11);
      if (x == 1)
        mesh.setAttribute(vertex.key(), 12);
    }
    Math::Matrix<Scalar> mass = Math::Matrix<Scalar>::Zero(3, 3);
    Math::Vector<Scalar> moments(3);
    const Real entries[3][3] = {{Real(1) / 6, Real(1) / 12, 0},
      {Real(1) / 12, Real(1) / 3, Real(1) / 12}, {0, Real(1) / 12, Real(1) / 6}};
    const Real nodalMoments[3] = {Real(1) / 4, Real(1) / 2, Real(1) / 4};
    for (size_t i = 0; i < 3; ++i)
    {
      moments[dofs[i]] = nodalMoments[i];
      for (size_t j = 0; j < 3; ++j)
        mass(dofs[i], dofs[j]) = entries[i][j];
    }
    PETSc::Variational::TrialFunction u(space);
    PETSc::Variational::TestFunction v(space);
    if constexpr (!Mixed)
    {
      Problem problem(u, v);
      problem = Integral(u, v) - Integral(known(FieldScalar(2) * phase), v) +
        DirichletBC(u, known(FieldScalar(2) * phase)).on(11);
      problem.assemble();
      const Math::Vector<Scalar> exact =
        Math::Vector<Scalar>::Constant(3, Scalar(2) * phase);
      checkOwnership(problem, mass, Scalar(2) * phase * moments, exact, {dofs[0]},
        [&] { EXPECT_LT((vector(u.getSolution().getData()) - exact).norm(), 1e-10); });
    }
    else
    {
      PETSc::Variational::TrialFunction p(space);
      PETSc::Variational::TestFunction q(space);
      Problem problem(u, v, p, q);
      problem = Integral(FieldScalar(2) * u, v) + Integral(p, v) + Integral(u, q) +
        Integral(FieldScalar(3) * p, q) - Integral(known(FieldScalar(0)), v) -
        Integral(known(FieldScalar(-5) * phase), q) +
        DirichletBC(u, known(phase)).on(11) +
        DirichletBC(p, known(FieldScalar(-2) * phase)).on(12);
      problem.assemble();
      const auto& offsets = problem.getTrialOffsets();
      Math::Matrix<Scalar> raw = Math::Matrix<Scalar>::Zero(6, 6);
      Math::Vector<Scalar> load = Math::Vector<Scalar>::Zero(6);
      Math::Vector<Scalar> exact(6);
      for (size_t i = 0; i < 2; ++i)
      {
        exact.segment(offsets[i], 3).setConstant(i == 0 ? phase : Scalar(-2) * phase);
        for (size_t j = 0; j < 2; ++j)
          raw.block(offsets[i], offsets[j], 3, 3) = Scalar(i == j ? i + 2 : 1) * mass;
      }
      load.segment(offsets[1], 3) = Scalar(-5) * phase * moments;
      checkOwnership(
        problem, raw, load, exact, {offsets[0] + dofs[0], offsets[1] + dofs[2]}, [&] {
          EXPECT_LT(
            (vector(u.getSolution().getData()) - exact.segment(offsets[0], 3)).norm(),
            1e-10);
          EXPECT_LT(
            (vector(p.getSolution().getData()) - exact.segment(offsets[1], 3)).norm(),
            1e-10);
        });
    }
  }
}

TEST(PETSc_ProblemLifecycle, SingleField)
{
  checkProblemLifecycle<false>();
}

TEST(PETSc_ProblemLifecycle, MixedFields)
{
  checkProblemLifecycle<true>();
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  ::testing::InitGoogleTest(&argc, argv);
  const int result = RUN_ALL_TESTS();
  PetscFinalize();
  return result;
}
