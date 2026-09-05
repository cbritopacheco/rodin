/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIFFUSIONFORM_H
#define RODIN_VARIATIONAL_DIFFUSIONFORM_H

#include <functional>
#include <vector>

#include "Rodin/Assembly/Default.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

#include "BilinearForm.h"

namespace Rodin::FormLanguage
{
  template <class Solution, class TrialFES, class TestFES, class Operator>
  struct Traits<Variational::DiffusionForm<Solution, TrialFES, TestFES, Operator>>
  {
    using SolutionType = Solution;
    using TrialFESType = TrialFES;
    using TestFESType = TestFES;
    using OperatorType = Operator;
  };
}

namespace Rodin::Variational
{
  /**
   * @brief Discrete diffusion form
   * @f$ a(u,v) = \int_\Omega \nabla u : \nabla v\,dx @f$.
   */
  template <class Solution, class TrialFES, class TestFES, class Operator>
  class DiffusionForm final : public BilinearFormBase<Operator>
  {
    using TrialMeshType = typename FormLanguage::Traits<TrialFES>::MeshType;
    using TestMeshType = typename FormLanguage::Traits<TestFES>::MeshType;
    using TrialContextType = typename FormLanguage::Traits<TrialMeshType>::ContextType;
    using TestContextType = typename FormLanguage::Traits<TestMeshType>::ContextType;

    public:
      using SolutionType = Solution;
      using TrialFESType = TrialFES;
      using TestFESType = TestFES;
      using OperatorType = Operator;
      using ScalarType = typename FormLanguage::Traits<OperatorType>::ScalarType;
      /**
       * @brief Cell kernel of the diffusion form.
       *
       * Evaluates the reference gradients once per polytope geometry and
       * reuses them across cells, pulling them back through the metric
       * @f$ J^{-1} J^{-T} @f$ at each cell.
       */
      class Kernel
      {
        public:
          using MatrixType = Math::Matrix<ScalarType>;

          Kernel(const TrialFES& trialFES, const TestFES& testFES)
            : m_trialFES(trialFES), m_testFES(testFES)
          {}

          void compute(MatrixType& out, const Geometry::Polytope& polytope) const
          {
            const size_t d = polytope.getDimension();
            const Index i = polytope.getIndex();
            const auto& trialFE = m_trialFES.get().getFiniteElement(d, i);
            const auto& testFE = m_testFES.get().getFiniteElement(d, i);
            const size_t order =
              trialFE.getOrder() == 0 || testFE.getOrder() == 0
                ? 0
                : trialFE.getOrder() + testFE.getOrder() - 2;
            const auto& qf = QF::PolytopeQuadratureFormula::get(
              order, polytope.getGeometry());
            const auto& quadrature = polytope.getQuadrature(qf);

            const bool rebuild =
              !m_cached || m_geometry != polytope.getGeometry()
              || m_order != order
              || m_trialCount != trialFE.getCount()
              || m_testCount != testFE.getCount();
            if (rebuild)
            {
              m_cached = true;
              m_geometry = polytope.getGeometry();
              m_order = order;
              m_trialCount = trialFE.getCount();
              m_testCount = testFE.getCount();
              m_trialGradients.assign(qf.getSize(), {});
              m_testGradients.assign(qf.getSize(), {});
              for (size_t qp = 0; qp < qf.getSize(); ++qp)
              {
                const auto& reference = qf.getPoint(qp);
                auto& trialGradients = m_trialGradients[qp];
                auto& testGradients = m_testGradients[qp];
                trialGradients.reserve(trialFE.getCount());
                testGradients.reserve(testFE.getCount());
                for (size_t tr = 0; tr < trialFE.getCount(); ++tr)
                  trialGradients.emplace_back(
                    trialFE.getBasis(tr).getGradient()(reference));
                for (size_t te = 0; te < testFE.getCount(); ++te)
                  testGradients.emplace_back(
                    testFE.getBasis(te).getGradient()(reference));
              }
            }

            out.resize(testFE.getCount(), trialFE.getCount());
            out.setZero();
            // The pullback of the physical gradients is folded into the metric
            // G = J^{-1} J^{-T}, so that only the trial gradients are mapped.
            m_scratch.resize(trialFE.getCount());
            for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
            {
              const auto& point = quadrature.getPoint(qp);
              const ScalarType weight = static_cast<ScalarType>(
                qf.getWeight(qp) * point.getDistortion());
              const auto& jacobianInverse = point.getJacobianInverse();
              const auto metric = jacobianInverse * jacobianInverse.transpose();
              const auto& trialGradients = m_trialGradients[qp];
              const auto& testGradients = m_testGradients[qp];
              for (size_t tr = 0; tr < trialFE.getCount(); ++tr)
                m_scratch[tr] = metric * trialGradients[tr];
              for (size_t te = 0; te < testFE.getCount(); ++te)
              {
                for (size_t tr = 0; tr < trialFE.getCount(); ++tr)
                  out(te, tr) +=
                    weight * Math::dot(testGradients[te], m_scratch[tr]);
              }
            }
          }

        private:
          std::reference_wrapper<const TrialFES> m_trialFES;
          std::reference_wrapper<const TestFES> m_testFES;
          mutable bool m_cached = false;
          mutable Geometry::Polytope::Type m_geometry =
            Geometry::Polytope::Type::Point;
          mutable size_t m_order = 0;
          mutable size_t m_trialCount = 0;
          mutable size_t m_testCount = 0;
          mutable std::vector<std::vector<Math::SpatialVector<ScalarType>>>
            m_trialGradients;
          mutable std::vector<std::vector<Math::SpatialVector<ScalarType>>>
            m_testGradients;
          mutable std::vector<Math::SpatialVector<ScalarType>> m_scratch;
      };

      using KernelType = Kernel;
      using Parent = BilinearFormBase<OperatorType>;
      using AssemblyType =
        typename Assembly::Default<TrialContextType, TestContextType>
          ::template Type<OperatorType, DiffusionForm>;

      DiffusionForm(
          const TrialFunction<SolutionType, TrialFESType>& u,
          const TestFunction<TestFESType>& v)
        : m_u(u), m_v(v)
      {
        assemble();
      }

      DiffusionForm(const DiffusionForm&) = default;
      DiffusionForm(DiffusionForm&&) = default;
      DiffusionForm& operator=(const DiffusionForm&) = default;
      DiffusionForm& operator=(DiffusionForm&&) = default;

      OperatorType& getOperator() override { return m_operator; }
      const OperatorType& getOperator() const override { return m_operator; }

      void assemble() override
      {
        m_assembly.execute(m_operator, *this);
      }

      const TrialFunction<SolutionType, TrialFESType>& getTrialFunction() const override
      {
        return m_u.get();
      }

      const TestFunction<TestFESType>& getTestFunction() const override
      {
        return m_v.get();
      }

      DiffusionForm* copy() const noexcept override { return new DiffusionForm(*this); }

    private:
      std::reference_wrapper<const TrialFunction<SolutionType, TrialFESType>> m_u;
      std::reference_wrapper<const TestFunction<TestFESType>> m_v;
      OperatorType m_operator;
      AssemblyType m_assembly;
  };

  template <class Solution, class TrialFES, class TestFES>
  DiffusionForm(const TrialFunction<Solution, TrialFES>&, const TestFunction<TestFES>&)
    -> DiffusionForm<
        Solution, TrialFES, TestFES,
        Math::SparseMatrix<
          typename FormLanguage::Mult<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>;
}

#endif
