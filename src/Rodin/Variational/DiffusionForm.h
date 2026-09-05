/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIFFUSIONFORM_H
#define RODIN_VARIATIONAL_DIFFUSIONFORM_H

#include <functional>

#include "Rodin/Assembly/Default.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Math/SparseMatrix.h"

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
