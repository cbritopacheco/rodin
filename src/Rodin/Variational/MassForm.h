/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MASSFORM_H
#define RODIN_VARIATIONAL_MASSFORM_H

#include <functional>

#include "Rodin/Assembly/Default.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Math/SparseMatrix.h"

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
