/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef RODIN_TESTS_CONVERGENCE_PETSC_LINEARELASTICITYREFINEMENT_H
#define RODIN_TESTS_CONVERGENCE_PETSC_LINEARELASTICITYREFINEMENT_H

#include <utility>
#include <type_traits>

#include "FieldConvergence.h"
#include "PETScLinearElasticityProblem.h"

namespace Rodin::Tests::Convergence
{
  /** @brief Vector verification along degree and combined paths.
   * @par Architecture
   * A mesh factory supplies local or distributed grids. Each measurement owns
   * a fresh vector workload. FieldConvergence checks every field and every
   * interval; polynomial controls have a separate absolute-error contract.
   */
  template <class MeshFactory>
  class PETScLinearElasticityRefinement
  {
    private:
      using Data = LinearElasticity::ManufacturedSolution;
      static constexpr size_t AssemblyOrder = 16;
      static constexpr size_t NormOrder = AssemblyOrder + 2;
      static constexpr Real SolverTolerance = 1e-13;
      static constexpr Real PatchTolerance = 1e-9;
      static constexpr Real SensitivityTolerance = 1e-6;
      // Finite-path acceptance policies shared with the native vector studies.
      static constexpr Real DegreeFloor = 0.25;
      static constexpr Real EffectiveL2Floor = 1.9, EffectiveH1Floor = 0.9;
      // Separated errors of the resolved patch solved with the omitted volumetric term.
      static constexpr Real OmittedVolumetricL2Floor = 1e-3;
      static constexpr Real OmittedVolumetricH1Floor = 1e-2;

    public:
      explicit PETScLinearElasticityRefinement(MeshFactory factory)
        : m_factory(std::move(factory))
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
        expectExact(solve<K>(mesh, Data::Field::Quadratic));
      }

      template <size_t K>
      void checkVolumetricControl() const
      {
        const auto mesh = m_factory(3);
        const auto field = Data::Field::Quadratic;
        expectExact(solve<K>(mesh, field));
        const auto error = solve<K>(mesh, field, true);
        ASSERT_TRUE(error.isFinite());
        EXPECT_GT(error.getL2(), OmittedVolumetricL2Floor);
        EXPECT_GT(error.getH1Seminorm(), OmittedVolumetricH1Floor);
      }

      void checkDegreeThreshold() const
      {
        const auto mesh = m_factory(2);
        const auto error = solve<1>(mesh, Data::Field::Quadratic);
        ASSERT_TRUE(error.isFinite());
        EXPECT_GT(error.getL2(), OmittedVolumetricL2Floor);
        EXPECT_GT(error.getH1Seminorm(), OmittedVolumetricH1Floor);
        expectExact(solve<2>(mesh, Data::Field::Quadratic));
      }

      template <size_t K>
      void checkAsymmetricTensorPatch() const
      {
        const auto mesh = m_factory(3);
        PETScLinearElasticityProblem<K, std::remove_cv_t<decltype(mesh)>> problem(
          mesh, Data::Field::AsymmetricAffine, AssemblyOrder);
        const auto error = problem.solve(
          false, SolverTolerance, NormOrder, [&](const auto& solution, const Data& data) {
            const auto jacobian = Variational::Jacobian(solution);
            const auto tensor = [&](
                                  const Variational::IntegrationPoint& ip, bool stress) {
              const auto gradient = jacobian(ip);
              auto value = gradient;
              Real trace = 0;
              for (size_t i = 0; i < gradient.rows(); ++i)
                trace += gradient(i, i);
              for (size_t i = 0; i < gradient.rows(); ++i)
                for (size_t j = 0; j < gradient.cols(); ++j)
                  value(i, j) = stress ? data.getLambda() * trace * (i == j) +
                      data.getMu() * (gradient(i, j) + gradient(j, i))
                                       : (gradient(i, j) + gradient(j, i)) / 2;
              return value;
            };
            const auto strain = [&](const Variational::IntegrationPoint& ip) {
              return tensor(ip, false);
            };
            const auto stress = [&](const Variational::IntegrationPoint& ip) {
              return tensor(ip, true);
            };
            EXPECT_LT(
              ErrorNorm::computeL2(
                mesh, strain,
                [&data](const Geometry::Point& p) { return data.strain(p); }, NormOrder),
              PatchTolerance);
            EXPECT_LT(
              ErrorNorm::computeL2(
                mesh, stress,
                [&data](const Geometry::Point& p) { return data.stress(p); }, NormOrder),
              PatchTolerance);
          });
        expectExact(error);
      }

      template <size_t K>
      void checkSensitivity() const
      {
        const auto mesh = m_factory(3);
        const auto baseline = solve<K>(mesh);
        for (const auto& refined :
          {solve<K>(mesh, Data::Field::Exponential, false, AssemblyOrder + 2),
            solve<K>(mesh, Data::Field::Exponential, false, AssemblyOrder,
              SolverTolerance, NormOrder + 2),
            solve<K>(mesh, Data::Field::Exponential, false, AssemblyOrder,
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
      ErrorNorms solve(const MeshType& mesh, Data::Field field = Data::Field::Exponential,
        bool wrong = false, size_t assemblyOrder = AssemblyOrder,
        Real tolerance = SolverTolerance, size_t normOrder = NormOrder) const
      {
        return PETScLinearElasticityProblem<K, MeshType>(mesh, field, assemblyOrder)
          .solve(wrong, tolerance, normOrder);
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
