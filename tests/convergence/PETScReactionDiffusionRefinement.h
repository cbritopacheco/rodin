/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_REACTIONDIFFUSIONREFINEMENT_H
#define RODIN_TESTS_CONVERGENCE_PETSC_REACTIONDIFFUSIONREFINEMENT_H

#include <utility>

#include "FieldConvergence.h"
#include "PETScReactionDiffusionProblem.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Componentwise verification along degree and combined paths.
   * @par Architecture
   * A mesh factory supplies local or distributed grids. Each measurement owns
   * a fresh two-field workload. FieldConvergence checks every field and every
   * interval; polynomial controls have a separate absolute-error contract.
   */
  template <class MeshFactory>
  class PETScReactionDiffusionRefinement
  {
    private:
      static constexpr size_t AssemblyOrder = 16;
      static constexpr size_t NormOrder = AssemblyOrder + 2;
      static constexpr Real SolverTolerance = 1e-13;
      static constexpr Real PatchTolerance = 1e-9;
      static constexpr Real SensitivityTolerance = 1e-6;
      // Finite-path acceptance policies shared with the native coupled study.
      static constexpr Real DegreeFloor = 0.25;
      static constexpr Real EffectiveL2Floor = 1.9, EffectiveH1Floor = 0.9;
      // Separated errors of the affine problem solved without its coupling.
      static constexpr Real WrongCouplingL2Floor = 1e-3;
      static constexpr Real WrongCouplingH1Floor = 1e-2;

    public:
      explicit PETScReactionDiffusionRefinement(MeshFactory factory)
        : m_factory(std::move(factory))
      {}

      void checkDegreePath() const
      {
        const auto mesh = m_factory(2);
        FieldConvergence<2> history;
        history.append(1, solve<1>(mesh));
        history.append(2, solve<2>(mesh));
        history.append(3, solve<3>(mesh));
        history.append(4, solve<4>(mesh));
        ASSERT_EQ(history.getSize(), 4u);
        history.expectExponentialFloor({DegreeFloor, DegreeFloor});
      }

      void checkCombinedPath() const
      {
        const auto coarse = m_factory(2), middle = m_factory(3), fine = m_factory(5);
        FieldConvergence<2> history;
        history.append(1, solve<1>(coarse));
        history.append(0.5, solve<2>(middle));
        history.append(0.25, solve<3>(fine));
        ASSERT_EQ(history.getSize(), 3u);
        history.expectAlgebraicFloor({EffectiveL2Floor, EffectiveH1Floor});
      }

      template <size_t K>
      void checkPatch() const
      {
        const auto mesh = m_factory(3);
        for (const auto& error : solve<K>(mesh, ReactionDiffusionData::Field::Quadratic))
          expectExact(error);
      }

      template <size_t K>
      void checkCouplingControl() const
      {
        const auto mesh = m_factory(3);
        for (const auto& error : solve<K>(mesh, ReactionDiffusionData::Field::Affine))
          expectExact(error);
        for (const auto& error :
          solve<K>(mesh, ReactionDiffusionData::Field::Affine, true))
        {
          ASSERT_TRUE(error.isFinite());
          EXPECT_GT(error.getL2(), WrongCouplingL2Floor);
          EXPECT_GT(error.getH1Seminorm(), WrongCouplingH1Floor);
        }
      }

      template <size_t K>
      void checkSensitivity() const
      {
        const auto mesh = m_factory(3);
        const auto baseline = solve<K>(mesh);
        for (const auto& refined :
          {solve<K>(mesh, ReactionDiffusionData::Field::Smooth, false, AssemblyOrder + 2),
            solve<K>(mesh, ReactionDiffusionData::Field::Smooth, false, AssemblyOrder,
              SolverTolerance, NormOrder + 2),
            solve<K>(mesh, ReactionDiffusionData::Field::Smooth, false, AssemblyOrder,
              SolverTolerance / 10)})
          for (size_t field = 0; field < 2; ++field)
          {
            SCOPED_TRACE(::testing::Message() << "field=" << field);
            ASSERT_TRUE(baseline[field].isFinite());
            ASSERT_TRUE(refined[field].isFinite());
            ASSERT_GT(baseline[field].getL2(), 0);
            ASSERT_GT(baseline[field].getH1Seminorm(), 0);
            EXPECT_LT(std::abs(refined[field].getL2() / baseline[field].getL2() - 1),
              SensitivityTolerance);
            EXPECT_LT(
              std::abs(
                refined[field].getH1Seminorm() / baseline[field].getH1Seminorm() - 1),
              SensitivityTolerance);
          }
      }

    private:
      template <size_t K, class MeshType>
      static std::array<ErrorNorms, 2> solve(const MeshType& mesh,
        ReactionDiffusionData::Field field = ReactionDiffusionData::Field::Smooth,
        bool omitCoupling = false, size_t assemblyOrder = AssemblyOrder,
        Real tolerance = SolverTolerance, size_t normOrder = NormOrder)
      {
        return PETScReactionDiffusionProblem<K, MeshType>(mesh, field, assemblyOrder)
          .solve(omitCoupling, tolerance, normOrder);
      }

      static void expectExact(const ErrorNorms& error)
      {
        ASSERT_TRUE(error.isFinite());
        EXPECT_LT(error.getL2(), PatchTolerance);
        EXPECT_LT(error.getH1Seminorm(), PatchTolerance);
      }

      MeshFactory m_factory;
  };
}

#endif
