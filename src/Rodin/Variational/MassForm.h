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
  /**
   * @brief Traits of a Variational::MassForm.
   *
   * @tparam Solution Solution variable type.
   * @tparam TrialFES Trial finite element space type.
   * @tparam TestFES Test finite element space type.
   * @tparam Operator Assembled operator type.
   */
  template <class Solution, class TrialFES, class TestFES, class Operator>
  struct Traits<Variational::MassForm<Solution, TrialFES, TestFES, Operator>>
  {
    /// @brief Solution variable type.
      using SolutionType = Solution;

    /// @brief Trial finite element space type.
      using TrialFESType = TrialFES;

    /// @brief Test finite element space type.
      using TestFESType = TestFES;

    /// @brief Assembled operator type.
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
      /// @brief Solution variable type.
      using SolutionType = Solution;

      /// @brief Trial finite element space type.
      using TrialFESType = TrialFES;

      /// @brief Test finite element space type.
      using TestFESType = TestFES;

      /// @brief Assembled operator type.
      using OperatorType = Operator;

      /// @brief Scalar value type of the operator.
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
          /// @brief Cell matrix type produced by the kernel.
          using MatrixType = Math::Matrix<ScalarType>;

          /**
           * @brief Constructs the kernel over the two spaces.
           * @param[in] trialFES Trial finite element space.
           * @param[in] testFES Test finite element space.
           */
          Kernel(const TrialFES& trialFES, const TestFES& testFES)
            : m_trialFES(trialFES),
              m_testFES(testFES)
          {}

          /**
           * @brief Computes the cell matrix of @p polytope.
           *
           * The reference bases are evaluated once per polytope geometry and
           * reused across cells.
           *
           * @param[in,out] out Cell matrix, resized to the local DOF counts.
           * @param[in] polytope Cell to integrate over.
           */
          void compute(MatrixType& out, const Geometry::Polytope& polytope) const
          {
            const size_t d = polytope.getDimension();
            const Index i = polytope.getIndex();
            const auto& trialFE = m_trialFES.get().getFiniteElement(d, i);
            const auto& testFE = m_testFES.get().getFiniteElement(d, i);
            const size_t order = trialFE.getOrder() + testFE.getOrder();
            const auto& qf =
              QF::PolytopeQuadratureFormula::get(order, polytope.getGeometry());
            const auto& quadrature = polytope.getQuadrature(qf);

            const bool rebuild = !m_cached || m_geometry != polytope.getGeometry() ||
              m_order != order || m_trialCount != trialFE.getCount() ||
              m_testCount != testFE.getCount();
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
                    matrix(te, tr) =
                      Math::dot(trialFE.getBasis(tr)(reference), testValue);
                }
              }
            }

            out.resize(testFE.getCount(), trialFE.getCount());
            out.setZero();
            for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
            {
              const auto& point = quadrature.getPoint(qp);
              const ScalarType weight =
                static_cast<ScalarType>(qf.getWeight(qp) * point.getDistortion());
              out += weight * m_reference[qp];
            }
          }

        private:
          std::reference_wrapper<const TrialFES> m_trialFES;
          std::reference_wrapper<const TestFES> m_testFES;
          mutable bool m_cached = false;
          mutable Geometry::Polytope::Type m_geometry = Geometry::Polytope::Type::Point;
          mutable size_t m_order = 0;
          mutable size_t m_trialCount = 0;
          mutable size_t m_testCount = 0;
          mutable std::vector<MatrixType> m_reference;
      };

      /// @brief Cell kernel of the form.
      using KernelType = Kernel;

      /// @brief Parent class.
      using Parent = BilinearFormBase<OperatorType>;

      /// @brief Assembly backend selected by the mesh contexts.
      using AssemblyType = typename Assembly::Default<TrialContextType,
        TestContextType>::template Type<OperatorType, MassForm>;

      /**
       * @brief Constructs the form over @p u and @p v and assembles it.
       * @param[in] u Trial function.
       * @param[in] v Test function.
       */
      MassForm(const TrialFunction<SolutionType, TrialFESType>& u,
        const TestFunction<TestFESType>& v)
        : m_u(u),
          m_v(v)
      {
        assemble();
      }

      /// @brief Copy constructor.
      MassForm(const MassForm&) = default;

      /// @brief Move constructor.
      MassForm(MassForm&&) = default;

      /// @brief Copy assignment.
      MassForm& operator=(const MassForm&) = default;

      /// @brief Move assignment.
      MassForm& operator=(MassForm&&) = default;

      OperatorType& getOperator() override
      {
        return m_operator;
      }
      const OperatorType& getOperator() const override
      {
        return m_operator;
      }

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

      MassForm* copy() const noexcept override
      {
        return new MassForm(*this);
      }

    private:
      std::reference_wrapper<const TrialFunction<SolutionType, TrialFESType>> m_u;
      std::reference_wrapper<const TestFunction<TestFESType>> m_v;
      OperatorType m_operator;
      AssemblyType m_assembly;
  };

  /// @brief Template argument deduction guide for MassForm.
  template <class Solution, class TrialFES, class TestFES>
  MassForm(const TrialFunction<Solution, TrialFES>&, const TestFunction<TestFES>&)
    -> MassForm<Solution, TrialFES, TestFES,
      Math::SparseMatrix<
        typename FormLanguage::Mult<typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>;
}

#endif
