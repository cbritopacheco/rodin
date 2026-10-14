/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Exact mass and boundary oracles for single-field PETSc targets. */
#ifndef RODIN_TESTS_PETSC_SINGLEFIELDTARGETASSEMBLY_H
#define RODIN_TESTS_PETSC_SINGLEFIELDTARGETASSEMBLY_H
#include <gtest/gtest.h>
#include <petsc.h>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>

namespace Rodin::Tests::Unit
{
  using namespace Geometry;
  using namespace Variational;

  /// @brief Gathers the owned entries of a PETSc vector without reading remote rows.
  inline Math::Vector<PetscScalar> gatherTargetVector(::Vec vector, MPI_Comm communicator)
  {
    PetscInt size = 0, first = 0, last = 0;
    EXPECT_EQ(VecGetSize(vector, &size), PETSC_SUCCESS);
    EXPECT_EQ(VecGetOwnershipRange(vector, &first, &last), PETSC_SUCCESS);
    Math::Vector<PetscScalar> values(size);
    values.setZero();
    std::vector<PetscInt> indices(static_cast<size_t>(last - first));
    std::iota(indices.begin(), indices.end(), first);
    if (first != last)
    {
      EXPECT_EQ(VecGetValues(vector, last - first, indices.data(), values.data() + first),
        PETSC_SUCCESS);
    }
    EXPECT_EQ(MPI_Allreduce(
                MPI_IN_PLACE, values.data(), size, MPIU_SCALAR, MPI_SUM, communicator),
      MPI_SUCCESS);
    return values;
  }

  /// @brief Labels nodal coordinates and components with separated scalar intervals.
  template <size_t Rank, class Scalar>
  auto singleFieldCoordinateFunction()
  {
    using Range = std::conditional_t<Rank == 0, Scalar,
      std::conditional_t<Rank == 1, Math::SpatialVector<Scalar>,
        Math::SpatialMatrix<Scalar>>>;
    return [](const Point& point) -> Range {
      const auto entry = [&](size_t component) -> Scalar {
        // Offset four keeps roundoff away from component-label interval boundaries.
        return Scalar(Real(16 * component) + Real(4) + point.x());
      };
      if constexpr (Rank == 0)
        return entry(0);
      else if constexpr (Rank == 1)
      {
        Range result(2);
        for (size_t component = 0; component < 2; ++component)
          result[component] = entry(component);
        return result;
      }
      else
      {
        Range result(2, 3);
        for (size_t row = 0; row < 2; ++row)
        {
          for (size_t col = 0; col < 3; ++col)
            result(row, col) = entry(3 * row + col);
        }
        return result;
      }
    };
  }

  /// @brief Creates distinct constant component values independently of the backend.
  template <size_t Rank, class Scalar>
  auto singleFieldBoundaryFunction(Scalar factor)
  {
    if constexpr (Rank == 0)
    {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return ComplexFunction(factor);
      else
        return RealFunction(factor);
    }
    else if constexpr (Rank == 1)
    {
      Math::SpatialVector<Scalar> constant(2);
      for (size_t c = 0; c < 2; ++c)
        constant[c] = Scalar(c + 1) * factor;
      return VectorFunction(size_t(2), [constant](const Point&) { return constant; });
    }
    else
    {
      Math::SpatialMatrix<Scalar> constant(2, 3);
      for (size_t row = 0; row < 2; ++row)
      {
        for (size_t col = 0; col < 3; ++col)
          constant(row, col) = Scalar(3 * row + col + 1) * factor;
      }
      return MatrixFunction(constant);
    }
  }

  /// @brief Exact nodal mass matrices distinguish every field component and target update.
  template <template <class, class> class Assembler, bool Quadratic, size_t Rank,
    class MeshType, class Scalar = PetscScalar>
  void checkSingleFieldTargets(MeshType& mesh, MPI_Comm communicator)
  {
    SCOPED_TRACE(::testing::Message() << "quadratic=" << Quadratic << ", rank=" << Rank);
    using Range = std::conditional_t<Rank == 0, Scalar,
      std::conditional_t<Rank == 1, Math::SpatialVector<Scalar>,
        Math::SpatialMatrix<Scalar>>>;
    constexpr size_t components = Rank == 0 ? 1 : Rank == 1 ? 2 : 6;
    constexpr size_t nodes = Quadratic ? 5 : 3;
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, -0.5);
      else
        return Scalar(1);
    }();
    const Scalar coefficient = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(2, 0.5);
      else
        return Scalar(2);
    }();
    auto space = [&] {
      if constexpr (Quadratic)
      {
        if constexpr (Rank == 0)
          return H1<2, Range, MeshType>(std::integral_constant<size_t, 2>{}, mesh);
        else if constexpr (Rank == 1)
          return H1<2, Range, MeshType>(std::integral_constant<size_t, 2>{}, mesh, 2);
        else
          return H1<2, Range, MeshType>(std::integral_constant<size_t, 2>{}, mesh, 2, 3);
      }
      else
      {
        if constexpr (Rank == 0)
          return P1<Range, MeshType>(mesh);
        else if constexpr (Rank == 1)
          return P1<Range, MeshType>(mesh, 2);
        else
          return P1<Range, MeshType>(mesh, 2, 3);
      }
    }();
    ASSERT_EQ(space.getSize(), components * nodes);
    PETSc::Variational::GridFunction coordinates(space);
    coordinates = singleFieldCoordinateFunction<Rank, Scalar>();
    coordinates.flush();
    const auto encoded = gatherTargetVector(coordinates.getData(), communicator);
    const Index absent = std::numeric_limits<Index>::max();
    std::vector<std::vector<Index>> dofs(components, std::vector<Index>(nodes, absent));
    for (Index dof = 0; dof < static_cast<Index>(encoded.size()); ++dof)
    {
      const Real code = PetscRealPart(encoded[dof]);
      EXPECT_LT(PetscAbsScalar(encoded[dof] - Scalar(code)), 1e-13);
      const size_t component = static_cast<size_t>(std::floor(code / Real(16)));
      ASSERT_LT(component, components);
      const Real x = code - Real(16 * component) - Real(4);
      const auto node = static_cast<size_t>(std::llround(x * Real(nodes - 1)));
      ASSERT_LT(node, nodes);
      ASSERT_NEAR(x, Real(node) / Real(nodes - 1), 1e-13);
      ASSERT_EQ(dofs[component][node], absent);
      dofs[component][node] = dof;
    }
    PETSc::Variational::TrialFunction u(space);
    PETSc::Variational::TestFunction v(space);
    const auto prescribed = singleFieldBoundaryFunction<Rank>(factor);
    PETSc::Variational::BilinearForm preassembled(u, v);
    preassembled = Integral(Scalar(0.5) * u, v);
    preassembled.assemble();
    PETSc::Variational::LinearForm preassembledLoad(v);
    preassembledLoad = Integral(Scalar(0.5) * prescribed, v);
    preassembledLoad.assemble();
    using System = PETSc::Math::LinearSystem;
    using ProblemType = Problem<System, decltype(u), decltype(v)>;
    using Body = typename ProblemType::ProblemBodyType;
    using Backend = Assembler<System, ProblemType>;
    Body body = Integral(coefficient * u, v) - Integral(coefficient * prescribed, v) +
      DirichletBC(u, prescribed).on(Attribute(11), Attribute(12));
    body.getBFs().add(preassembled);
    body.getLFs().add(preassembledLoad);
    Assembly::ProblemAssemblyInput<Body, decltype(u), decltype(v)> input(body, u, v);
    System system(communicator);
    std::unique_ptr<Backend> backend;
    {
      Backend original;
      Backend copied(original);
      Backend moved(std::move(copied));
      backend.reset(moved.copy());
    }
    const PetscInt size = static_cast<PetscInt>(space.getSize());
    std::vector<PetscInt> indices(static_cast<size_t>(size));
    std::iota(indices.begin(), indices.end(), PetscInt(0));
    auto matrix = [&] {
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> result(
        size, size);
      result.setZero();
      PetscInt first = 0, last = 0;
      EXPECT_EQ(MatGetOwnershipRange(system.getOperator(), &first, &last), PETSC_SUCCESS);
      if (first != last)
      {
        EXPECT_EQ(MatGetValues(system.getOperator(), last - first, indices.data() + first,
                    size, indices.data(), result.data() + first * size),
          PETSC_SUCCESS);
      }
      EXPECT_EQ(MPI_Allreduce(MPI_IN_PLACE, result.data(), size * size, MPIU_SCALAR,
                  MPI_SUM, communicator),
        MPI_SUCCESS);
      return Math::Matrix<Scalar>(result);
    };
    Math::Matrix<Real> nodal = Math::Matrix<Real>::Zero(nodes, nodes);
    // Element mass on a half-unit interval: P1 h/6 and P2 h/30 times these tables.
    constexpr Real linear[2][2] = {{2, 1}, {1, 2}};
    constexpr Real quadratic[3][3] = {{4, 2, -1}, {2, 16, 2}, {-1, 2, 4}};
    constexpr size_t degree = Quadratic ? 2 : 1;
    for (size_t cell = 0; cell < 2; ++cell)
    {
      for (size_t i = 0; i <= degree; ++i)
      {
        for (size_t j = 0; j <= degree; ++j)
        {
          nodal(degree * cell + i, degree * cell + j) +=
            Real(0.5) * (Quadratic ? quadratic[i][j] / Real(30) : linear[i][j] / Real(6));
        }
      }
    }
    const Real linearMoments[3] = {Real(1) / 4, Real(1) / 2, Real(1) / 4};
    const Real quadraticMoments[5] = {
      Real(1) / 12, Real(1) / 3, Real(1) / 6, Real(1) / 3, Real(1) / 12};
    Math::Matrix<Scalar> raw = Math::Matrix<Scalar>::Zero(size, size);
    Math::Vector<Scalar> exact(size), load(size);
    std::vector<Index> fixed;
    const Scalar totalCoefficient = coefficient + Scalar(0.5);
    for (size_t component = 0; component < components; ++component)
    {
      fixed.push_back(dofs[component].front());
      fixed.push_back(dofs[component].back());
      for (size_t i = 0; i < nodes; ++i)
      {
        const Index row = dofs[component][i];
        ASSERT_NE(row, absent);
        exact[row] = Scalar(component + 1) * factor;
        load[row] = totalCoefficient * exact[row] *
          Scalar(Quadratic ? quadraticMoments[i] : linearMoments[i]);
        for (size_t j = 0; j < nodes; ++j)
          raw(row, dofs[component][j]) = totalCoefficient * nodal(i, j);
      }
    }
    ASSERT_LT((raw * exact - load).norm(), 1e-12);
    auto full = raw, rowOnly = raw;
    auto fullLoad = load, rowLoad = load;
    for (const auto dof : fixed)
    {
      fullLoad -= raw.col(dof) * exact[dof];
      full.row(dof).setZero();
      full.col(dof).setZero();
      full(dof, dof) = 1;
      rowOnly.row(dof).setZero();
      rowOnly(dof, dof) = 1;
      rowLoad[dof] = exact[dof];
    }
    for (const auto dof : fixed)
      fullLoad[dof] = exact[dof];
    for (size_t repeat = 0; repeat < 2; ++repeat)
    {
      SCOPED_TRACE(repeat);
      backend->execute(system, input);
      EXPECT_LT((matrix() - full).norm(), 1e-11);
      EXPECT_LT(
        (gatherTargetVector(system.getVector(), communicator) - fullLoad).norm(), 1e-11);
      const Math::Vector<Scalar> recovered =
        matrix().fullPivLu().solve(gatherTargetVector(system.getVector(), communicator));
      EXPECT_LT((recovered - exact).norm(), 1e-10);
      ASSERT_EQ(VecSet(system.getVector(), Scalar(11) * factor), PETSC_SUCCESS);
      backend->execute(system, input, AssemblyTarget::LHS);
      EXPECT_LT((matrix() - rowOnly).norm(), 1e-11);
      EXPECT_LT((gatherTargetVector(system.getVector(), communicator) -
                  Math::Vector<Scalar>::Constant(size, Scalar(11) * factor))
                  .norm(),
        1e-11);
      ASSERT_EQ(MatScale(system.getOperator(), 3), PETSC_SUCCESS);
      backend->execute(system, input, AssemblyTarget::RHS);
      EXPECT_LT((matrix() - Scalar(3) * rowOnly).norm(), 1e-11);
      EXPECT_LT(
        (gatherTargetVector(system.getVector(), communicator) - rowLoad).norm(), 1e-11);
    }
  }

  /// @brief Instantiates scalar, vector and rectangular matrix ranges in both FE families.
  template <template <class, class> class Assembler, class MeshType>
  void checkAllSingleFieldTargets(MeshType& mesh, MPI_Comm communicator)
  {
    checkSingleFieldTargets<Assembler, false, 0>(mesh, communicator);
    checkSingleFieldTargets<Assembler, false, 1>(mesh, communicator);
    checkSingleFieldTargets<Assembler, false, 2>(mesh, communicator);
    checkSingleFieldTargets<Assembler, true, 0>(mesh, communicator);
    checkSingleFieldTargets<Assembler, true, 1>(mesh, communicator);
    checkSingleFieldTargets<Assembler, true, 2>(mesh, communicator);
  }
}
#endif
