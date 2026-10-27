/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_NONLINEAR_REFINEMENT_H
#define RODIN_TESTS_CONVERGENCE_PETSC_NONLINEAR_REFINEMENT_H

#include <type_traits>

#include "PETScNonlinearPoisson.h"
#include "FieldConvergence.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Shared degree/combined-refinement oracles for the SNES workload.
   * @par Architecture
   * The mesh factory supplies local or distributed grids. Every measurement
   * constructs a fresh workload, preserving its fixed-layout solver contract.
   * All adjacent intervals are checked; an effective hp slope is not a
   * fixed-degree approximation-order assertion.
   */
  template <class MeshFactory>
  class PETScNonlinearRefinement
  {
    private:
      // Dimensionless numerical policies shared by both refinement axes.
      static constexpr size_t AssemblyOrder = 16;
      static constexpr size_t NormOrder = AssemblyOrder + 2;
      static constexpr size_t RefinedAssemblyOrder = AssemblyOrder + 2;
      static constexpr size_t RefinedNormOrder = NormOrder + 2;
      static constexpr Real NonlinearTolerance = 1e-11;
      static constexpr Real RefinedNonlinearTolerance = NonlinearTolerance / 10;
      static constexpr Real SensitivityTolerance = 1e-6; // Relative contamination budget.
      static constexpr Real TangentTolerance = 1e-6; // Centered-difference defect budget.
      static constexpr Real WrongTangentFloor = 1e-3; // Separated incorrect derivative.
      static constexpr Real ControlAmplitude = 4; // Resolved incorrect-reaction oracle.
      static constexpr Real ControlL2 = 0.05, ControlH1 = 0.2;
      // Finite-resolution policies, not universal asymptotic constants.
      static constexpr Real DegreeReduction = 0.1;
      static constexpr Real EffectiveL2Floor = 1.9, EffectiveH1Floor = 0.9;

    public:
      explicit PETScNonlinearRefinement(MeshFactory factory)
        : m_factory(std::move(factory))
      {}

      void checkDegreePath() const
      {
        const auto mesh = m_factory(3);
        FieldConvergence<1> history;
        history.append(1, {solve<1>(mesh)});
        history.append(2, {solve<2>(mesh)});
        history.append(3, {solve<3>(mesh)});
        history.append(4, {solve<4>(mesh)});
        expectImprovement(history, false);
      }

      void checkCombinedPath() const
      {
        FieldConvergence<1> history;
        const auto mesh1 = m_factory(2), mesh2 = m_factory(3), mesh3 = m_factory(5);
        history.append(1, {solve<1>(mesh1)});
        history.append(0.5, {solve<2>(mesh2)});
        history.append(0.25, {solve<3>(mesh3)});
        expectImprovement(history, true);
      }

      template <size_t HighestDegree>
      void checkTangents() const
      {
        static_assert(HighestDegree == 3 || HighestDegree == 4);
        const auto mesh = m_factory(3);
        expectTangent<1>(mesh);
        expectTangent<2>(mesh);
        expectTangent<3>(mesh);
        if constexpr (HighestDegree == 4)
          expectTangent<4>(mesh);
        EXPECT_GT(
          (PETScNonlinearPoissonProblem<HighestDegree, std::remove_cv_t<decltype(mesh)>>(
            mesh, AssemblyOrder, 1, false, true)
              .tangentDefect()),
          WrongTangentFloor);
      }

      template <size_t K>
      void checkReactionControl(size_t n) const
      {
        const auto mesh = m_factory(n);
        const auto correct = solve<K>(mesh, AssemblyOrder, ControlAmplitude);
        const auto incorrect = solve<K>(mesh, AssemblyOrder, ControlAmplitude, true);
        EXPECT_LT(correct.getL2(), ControlL2);
        EXPECT_LT(correct.getH1Seminorm(), ControlH1);
        EXPECT_GT(incorrect.getL2(), ControlL2);
        EXPECT_GT(incorrect.getH1Seminorm(), ControlH1);
      }

      template <size_t K>
      void checkSensitivity() const
      {
        const auto mesh = m_factory(3);
        const auto baseline = solve<K>(mesh);
        // Separate assembly, norm quadrature and SNES stopping controls.
        for (const auto& refined :
          {solve<K>(mesh, RefinedAssemblyOrder),
            solve<K>(mesh, AssemblyOrder, 1, false, NonlinearTolerance, RefinedNormOrder),
            solve<K>(mesh, AssemblyOrder, 1, false, RefinedNonlinearTolerance)})
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
      static ErrorNorms solve(const MeshType& mesh, size_t assemblyOrder = AssemblyOrder,
        Real amplitude = 1, bool omitCubic = false, Real tolerance = NonlinearTolerance,
        size_t normOrder = NormOrder)
      {
        PETScNonlinearPoissonProblem<K, MeshType> problem(
          mesh, assemblyOrder, amplitude, omitCubic);
        return problem.solve(tolerance, normOrder, [](const auto&, const auto&) {});
      }

      template <size_t K, class MeshType>
      static void expectTangent(const MeshType& mesh)
      {
        EXPECT_LT((PETScNonlinearPoissonProblem<K, MeshType>(mesh, AssemblyOrder)
                      .tangentDefect()),
          TangentTolerance);
      }

      static void expectImprovement(const FieldConvergence<1>& history, bool combined)
      {
        ASSERT_EQ(history.getSize(), combined ? 3u : 4u);
        if (combined)
          history.expectAlgebraicFloor({EffectiveL2Floor, EffectiveH1Floor});
        else
          history.expectExponentialFloor({DegreeReduction, DegreeReduction});
      }

      MeshFactory m_factory;
  };
}

#endif
