/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Shared exact oracles for local and distributed mixed PETSc targets. */
#ifndef RODIN_TESTS_PETSC_MIXEDTARGETASSEMBLY_H
#define RODIN_TESTS_PETSC_MIXEDTARGETASSEMBLY_H
#include <gtest/gtest.h>
#include <petsc.h>
#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/PETSc.h>

namespace Rodin::Tests::Unit
{
  using namespace Geometry;
  using namespace Variational;

  /// @brief Mixed local, nonlocal and preassembled terms obey separate target updates.
  template <template <class, class> class Assembler, bool Affine = false, class Space,
    class Scalar = PetscScalar>
  void checkMixedTargetAssembly(Space& space,
    const std::vector<std::pair<Real, Index>>& positions, MPI_Comm communicator)
  {
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, -0.5);
      else
        return Scalar(1);
    }();
    ASSERT_EQ(space.getSize(), 3u);
    ASSERT_EQ(positions.size(), 3u);
    for (size_t i = 0; i < positions.size(); ++i)
      ASSERT_EQ(positions[i].first, Real(i) / Real(2));
    PETSc::Variational::TrialFunction u(space), p(space);
    PETSc::Variational::TestFunction v(space), q(space);
    using System = PETSc::Math::LinearSystem;
    using ProblemType =
      Problem<System, decltype(u), decltype(v), decltype(p), decltype(q)>;
    using Body = typename ProblemType::ProblemBodyType;
    using Backend = Assembler<System, ProblemType>;
    const auto kernel00 = [](const Point&, const Point&) { return Scalar(0.5); };
    const auto kernel01 = [](const Point&, const Point&) { return Scalar(-1); };
    const auto kernel10 = [](const Point&, const Point&) { return Scalar(2); };
    const auto kernel11 = [](const Point&, const Point&) { return Scalar(0.25); };
    PETSc::Variational::BilinearForm preassembled(u, q);
    preassembled = Integral(Scalar(1.25) * u, q);
    preassembled.assemble();
    PETSc::Variational::LinearForm preassembledLoad(q);
    preassembledLoad = Integral(Scalar(3) * factor * q);
    preassembledLoad.assemble();
    const auto boundaryValue = [factor](Scalar value) {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return ComplexFunction(Scalar(value) * factor);
      else
        return RealFunction(Scalar(value) * factor);
    };
    Body body = Integral(Scalar(2) * u, v) + Integral(Scalar(3) * p, v) - Integral(u, q) +
      Integral(Scalar(5) * p, q) + Integral(Potential(kernel00, u), v) +
      Integral(Potential(kernel01, p), v) + Integral(Potential(kernel10, u), q) +
      Integral(Potential(kernel11, p), q) - Integral(Scalar(-1.5) * factor * v) -
      Integral(Scalar(-11.25) * factor * q);
    const Scalar coefficient = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(0.25, 0.5);
      else
        return Scalar(0.25);
    }();
    if constexpr (Affine)
    {
      body.getDBCs().add(DirichletBC(
        u, coefficient * p, boundaryValue(Scalar(1) + Scalar(2) * coefficient))
                           .on(11));
    }
    else
      body.getDBCs().add(DirichletBC(u, boundaryValue(1)).on(11));
    body.getDBCs().add(DirichletBC(p, boundaryValue(-2)).on(12));
    body.getBFs().add(preassembled);
    body.getLFs().add(preassembledLoad);
    Tuple<std::reference_wrapper<decltype(u)>, std::reference_wrapper<decltype(p)>>
      trials(std::ref(u), std::ref(p));
    Tuple<std::reference_wrapper<decltype(v)>, std::reference_wrapper<decltype(q)>> tests(
      std::ref(v), std::ref(q));
    std::array<size_t, 2> trialOffsets{0, 3}, testOffsets{0, 3};
    boost::bimap<FormLanguage::Base::UUID, size_t> trialMap, testMap;
    trialMap.insert({u.getUUID(), 0});
    trialMap.insert({p.getUUID(), 1});
    testMap.insert({v.getUUID(), 0});
    testMap.insert({q.getUUID(), 1});
    Assembly::ProblemAssemblyInput<Body, decltype(u), decltype(v), decltype(p),
      decltype(q)>
      input(body, trials, tests, trialOffsets, testOffsets, trialMap, testMap, 6, 6);
    System system(communicator);
    std::unique_ptr<Backend> backend;
    {
      Backend original;
      Backend copied(original);
      Backend moved(std::move(copied));
      backend.reset(moved.copy());
    }
    constexpr PetscInt size = 6;
    std::array<PetscInt, size> indices;
    std::iota(indices.begin(), indices.end(), PetscInt(0));
    auto matrix = [&] {
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> values(
        size, size);
      values.setZero();
      PetscInt first = 0, last = 0;
      EXPECT_EQ(MatGetOwnershipRange(system.getOperator(), &first, &last), PETSC_SUCCESS);
      if (first != last)
      {
        EXPECT_EQ(MatGetValues(system.getOperator(), last - first, indices.data() + first,
                    size, indices.data(), values.data() + first * size),
          PETSC_SUCCESS);
      }
      EXPECT_EQ(MPI_Allreduce(MPI_IN_PLACE, values.data(), size * size, MPIU_SCALAR,
                  MPI_SUM, communicator),
        MPI_SUCCESS);
      return Math::Matrix<Scalar>(values);
    };
    auto vector = [&] {
      Math::Vector<Scalar> values(size);
      values.setZero();
      PetscInt first = 0, last = 0;
      EXPECT_EQ(VecGetOwnershipRange(system.getVector(), &first, &last), PETSC_SUCCESS);
      if (first != last)
      {
        EXPECT_EQ(VecGetValues(system.getVector(), last - first, indices.data() + first,
                    values.data() + first),
          PETSC_SUCCESS);
      }
      EXPECT_EQ(MPI_Allreduce(
                  MPI_IN_PLACE, values.data(), size, MPIU_SCALAR, MPI_SUM, communicator),
        MPI_SUCCESS);
      return values;
    };
    Math::Matrix<Scalar> raw(size, size);
    Math::Vector<Scalar> load(size), exact(size);
    const Real mass[3][3] = {{Real(1) / 6, Real(1) / 12, 0},
      {Real(1) / 12, Real(1) / 3, Real(1) / 12}, {0, Real(1) / 12, Real(1) / 6}};
    const Real moments[3] = {Real(1) / 4, Real(1) / 2, Real(1) / 4};
    const Real local[2][2] = {{2, 3}, {Real(0.25), 5}};
    const Real nonlocal[2][2] = {{Real(0.5), -1}, {2, Real(0.25)}};
    for (size_t testBlock = 0; testBlock < 2; ++testBlock)
    {
      for (size_t i = 0; i < positions.size(); ++i)
      {
        const Index row = 3 * testBlock + positions[i].second;
        load[row] =
          Scalar(moments[i] * (testBlock == 0 ? Real(-1.5) : Real(-8.25))) * factor;
        exact[row] = Scalar(testBlock == 0 ? 1 : -2) * factor;
        for (size_t trialBlock = 0; trialBlock < 2; ++trialBlock)
        {
          for (size_t j = 0; j < positions.size(); ++j)
          {
            const Index col = 3 * trialBlock + positions[j].second;
            raw(row, col) = Scalar(local[testBlock][trialBlock] * mass[i][j] +
              nonlocal[testBlock][trialBlock] * moments[i] * moments[j]);
          }
        }
      }
    }
    ASSERT_LT((raw * exact - load).norm(), 1e-13);
    const std::array<Index, 2> fixed{
      positions.front().second, 3 + positions.back().second};
    auto full = raw;
    auto rowOnly = raw;
    auto fullLoad = load;
    auto rowLoad = load;
    if constexpr (Affine)
    {
      const Index slave = fixed[0], master = 3 + positions.front().second;
      Math::Matrix<Scalar> projection = Math::Matrix<Scalar>::Identity(size, size);
      projection.row(slave).setZero();
      projection(slave, master) = coefficient;
      projection.row(fixed[1]).setZero();
      Math::Vector<Scalar> lift = Math::Vector<Scalar>::Zero(size);
      lift[slave] = (Scalar(1) + Scalar(2) * coefficient) * factor;
      lift[fixed[1]] = exact[fixed[1]];
      full = projection.adjoint() * raw * projection;
      fullLoad = projection.adjoint() * (load - raw * lift);
      full(slave, slave) = 1;
      full(slave, master) = -coefficient;
      fullLoad[slave] = lift[slave];
      full(fixed[1], fixed[1]) = 1;
      fullLoad[fixed[1]] = exact[fixed[1]];
    }
    else
    {
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
    }
    for (size_t repeat = 0; repeat < 2; ++repeat)
    {
      SCOPED_TRACE(repeat);
      backend->execute(system, input);
      EXPECT_LT((matrix() - full).norm(), 1e-12);
      EXPECT_LT((vector() - fullLoad).norm(), 1e-12);
      const Math::Vector<Scalar> recovered = matrix().fullPivLu().solve(vector());
      EXPECT_LT((recovered - exact).norm(), 1e-11);
      if constexpr (Affine)
      {
        for (const auto target : {AssemblyTarget::LHS, AssemblyTarget::RHS})
        {
          EXPECT_THROW(backend->execute(system, input, target), Alert::Exception);
          EXPECT_LT((matrix() - full).norm(), 1e-12);
          EXPECT_LT((vector() - fullLoad).norm(), 1e-12);
        }
        continue;
      }
      ASSERT_EQ(VecSet(system.getVector(), Scalar(11) * factor), PETSC_SUCCESS);
      backend->execute(system, input, AssemblyTarget::LHS);
      EXPECT_LT((matrix() - rowOnly).norm(), 1e-12);
      EXPECT_LT(
        (vector() - Math::Vector<Scalar>::Constant(size, Scalar(11) * factor)).norm(),
        1e-12);
      ASSERT_EQ(MatScale(system.getOperator(), 3), PETSC_SUCCESS);
      backend->execute(system, input, AssemblyTarget::RHS);
      EXPECT_LT((matrix() - Scalar(3) * rowOnly).norm(), 1e-12);
      EXPECT_LT((vector() - rowLoad).norm(), 1e-12);
    }
  }
}

#endif
