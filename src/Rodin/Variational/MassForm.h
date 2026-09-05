/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MASSFORM_H
#define RODIN_VARIATIONAL_MASSFORM_H

#include <functional>
#include <vector>

#include "Rodin/Assembly/Default.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

#include "BilinearForm.h"

namespace Rodin::FormLanguage
{
  template <class Solution, class TrialFES, class TestFES, class Operator>
  struct Traits<Variational::MassForm<Solution, TrialFES, TestFES, Operator>>
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
   * @brief Discrete mass form @f$ m(u,v) = \int_\Omega u \cdot v\,dx @f$.
   */
  template <class Solution, class TrialFES, class TestFES, class Operator>
  class MassForm final : public BilinearFormBase<Operator>
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
       * @brief Cell kernel of the mass form.
       *
       * Evaluates the reference bases once per polytope geometry and reuses
       * them across cells; only the quadrature weights and the distortion of
       * the map change from one cell to the next.
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
            const size_t order = trialFE.getOrder() + testFE.getOrder();
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
              m_reference.resize(qf.getSize());
              for (size_t qp = 0; qp < qf.getSize(); ++qp)
              {
                auto& matrix = m_reference[qp];
                matrix.resize(testFE.getCount(), trialFE.getCount());
                const auto& reference = qf.getPoint(qp);
                for (size_t te = 0; te < testFE.getCount(); ++te)
                {
                  const auto testValue = testFE.getBasis(te)(reference);
                  for (size_t tr = 0; tr < trialFE.getCount(); ++tr)
                    matrix(te, tr) = Math::dot(
                      trialFE.getBasis(tr)(reference), testValue);
                }
              }
            }

            out.resize(testFE.getCount(), trialFE.getCount());
            out.setZero();
            for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
            {
              const auto& point = quadrature.getPoint(qp);
              const ScalarType weight = static_cast<ScalarType>(
                qf.getWeight(qp) * point.getDistortion());
              out += weight * m_reference[qp];
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
          mutable std::vector<MatrixType> m_reference;
      };

      using KernelType = Kernel;
      using Parent = BilinearFormBase<OperatorType>;
      using AssemblyType =
        typename Assembly::Default<TrialContextType, TestContextType>
          ::template Type<OperatorType, MassForm>;

      MassForm(
          const TrialFunction<SolutionType, TrialFESType>& u,
          const TestFunction<TestFESType>& v)
        : m_u(u), m_v(v)
      {
        assemble();
      }

      MassForm(const MassForm&) = default;
      MassForm(MassForm&&) = default;
      MassForm& operator=(const MassForm&) = default;
      MassForm& operator=(MassForm&&) = default;

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

      MassForm* copy() const noexcept override { return new MassForm(*this); }

    private:
      std::reference_wrapper<const TrialFunction<SolutionType, TrialFESType>> m_u;
      std::reference_wrapper<const TestFunction<TestFESType>> m_v;
      OperatorType m_operator;
      AssemblyType m_assembly;
  };

  template <class Solution, class TrialFES, class TestFES>
  MassForm(const TrialFunction<Solution, TrialFES>&, const TestFunction<TestFES>&)
    -> MassForm<
        Solution, TrialFES, TestFES,
        Math::SparseMatrix<
          typename FormLanguage::Mult<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>;
}

#endif
