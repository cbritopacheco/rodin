/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_DIFFUSIONREFINEMENT_H
#define RODIN_TESTS_CONVERGENCE_PETSC_DIFFUSIONREFINEMENT_H

#include <utility>

#include "FieldConvergence.h"
#include "PETScDiffusionProblem.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Scalar verification along degree and combined paths.
   * @par Architecture
   * A mesh factory supplies local or distributed grids. Each measurement owns
   * a fresh scalar workload. FieldConvergence checks every field and every
   * interval; polynomial controls have a separate absolute-error contract.
   */
  template <class MeshFactory>
  class PETScDiffusionRefinement
  {
    private:
      static constexpr size_t AssemblyOrder = 16;
      static constexpr size_t NormOrder = AssemblyOrder + 2;
      static constexpr Real SolverTolerance = 1e-13;
      static constexpr Real PatchTolerance = 1e-9;
      static constexpr Real SensitivityTolerance = 1e-6;
      // Finite-path acceptance policies shared with the native scalar studies.
      static constexpr Real DegreeFloor = 0.25;
      static constexpr Real EffectiveL2Floor = 1.9, EffectiveH1Floor = 0.9;
      // Separated errors of the resolved patch solved with the wrong coefficient.
      static constexpr Real WrongOperatorL2Floor = 1e-3;
      static constexpr Real WrongOperatorH1Floor = 1e-2;

    public:
      PETScDiffusionRefinement(MeshFactory factory, bool poisson)
        : m_factory(std::move(factory)),
          m_poisson(poisson)
      {}

      void checkDegreePath() const
      {
        const auto mesh = m_factory(2);
        FieldConvergence<1> history;
        history.append(1, {solve<1>(mesh)});
        history.append(2, {solve<2>(mesh)});
        history.append(3, {solve<3>(mesh)});
        history.append(4, {solve<4>(mesh)});
        ASSERT_EQ(history.getSize(), 4u);
        history.expectExponentialFloor({DegreeFloor, DegreeFloor});
      }

      void checkCombinedPath() const
      {
        const auto coarse = m_factory(2), middle = m_factory(3), fine = m_factory(5);
        FieldConvergence<1> history;
        history.append(1, {solve<1>(coarse)});
        history.append(0.5, {solve<2>(middle)});
        history.append(0.25, {solve<3>(fine)});
        ASSERT_EQ(history.getSize(), 3u);
        history.expectAlgebraicFloor({EffectiveL2Floor, EffectiveH1Floor});
      }

      template <size_t K>
      void checkPatch() const
      {
        const auto mesh = m_factory(3);
        expectExact(solve<K>(mesh, ConductivityData::Field::Quadratic));
      }

      template <size_t K>
      void checkOperatorControl() const
      {
        const auto mesh = m_factory(3);
        const auto field = m_poisson ? ConductivityData::Field::Quadratic
                                     : ConductivityData::Field::Affine;
        expectExact(solve<K>(mesh, field));
        const auto error = solve<K>(mesh, field, true);
        ASSERT_TRUE(error.isFinite());
        EXPECT_GT(error.getL2(), WrongOperatorL2Floor);
        EXPECT_GT(error.getH1Seminorm(), WrongOperatorH1Floor);
      }

      void checkDegreeThreshold() const
      {
        const auto mesh = m_factory(2);
        const auto error = solve<1>(mesh, ConductivityData::Field::Quadratic);
        ASSERT_TRUE(error.isFinite());
        EXPECT_GT(error.getL2(), WrongOperatorL2Floor);
        EXPECT_GT(error.getH1Seminorm(), WrongOperatorH1Floor);
        expectExact(solve<2>(mesh, ConductivityData::Field::Quadratic));
      }

      template <size_t K>
      void checkSensitivity() const
      {
        const auto mesh = m_factory(3);
        const auto baseline = solve<K>(mesh);
        for (const auto& refined :
          {solve<K>(mesh, ConductivityData::Field::Exponential, false, AssemblyOrder + 2),
            solve<K>(mesh, ConductivityData::Field::Exponential, false, AssemblyOrder,
              SolverTolerance, NormOrder + 2),
            solve<K>(mesh, ConductivityData::Field::Exponential, false, AssemblyOrder,
              SolverTolerance / 10)})
        {
          ASSERT_TRUE(baseline.isFinite());
          ASSERT_TRUE(refined.isFinite());
          ASSERT_GT(baseline.getL2(), 0);
          ASSERT_GT(baseline.getH1Seminorm(), 0);
          EXPECT_LT(
            std::abs(refined.getL2() / baseline.getL2() - 1), SensitivityTolerance);
          EXPECT_LT(std::abs(refined.getH1Seminorm() / baseline.getH1Seminorm() - 1),
            SensitivityTolerance);
        }
      }

    private:
      template <size_t K, class MeshType>
      ErrorNorms solve(const MeshType& mesh,
        ConductivityData::Field field = ConductivityData::Field::Exponential,
        bool wrong = false, size_t assemblyOrder = AssemblyOrder,
        Real tolerance = SolverTolerance, size_t normOrder = NormOrder) const
      {
        return PETScDiffusionProblem<K, MeshType>(mesh, m_poisson, field, assemblyOrder)
          .solve(wrong, tolerance, normOrder);
      }

      static void expectExact(const ErrorNorms& error)
      {
        ASSERT_TRUE(error.isFinite());
        EXPECT_LT(error.getL2(), PatchTolerance);
        EXPECT_LT(error.getH1Seminorm(), PatchTolerance);
      }

      MeshFactory m_factory;
      bool m_poisson;
  };
}

#endif
