/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MASSFORM_H
#define RODIN_VARIATIONAL_MASSFORM_H

#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

#include "Rodin/Assembly/Default.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/Region.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

#include "BilinearForm.h"
#include "Function.h"
#include "IntegrationPoint.h"
#include "RealFunction.h"

namespace Rodin::FormLanguage
{
  /// @brief Traits of a Variational::MassForm.
  template <class Solution, class TrialFES, class TestFES, class Operator,
    class CoefficientDerived>
  struct Traits<
    Variational::MassForm<Solution, TrialFES, TestFES, Operator, CoefficientDerived>>
  {
      /// @brief Solution type of the trial function.
      using SolutionType = Solution;

      /// @brief Trial finite element space type.
      using TrialFESType = TrialFES;

      /// @brief Test finite element space type.
      using TestFESType = TestFES;

      /// @brief Assembled operator type.
      using OperatorType = Operator;
  };

  /// @brief Marks Variational::MassForm as a named form.
  template <class Solution, class TrialFES, class TestFES, class Operator,
    class CoefficientDerived>
  struct IsNamedForm<
    Variational::MassForm<Solution, TrialFES, TestFES, Operator, CoefficientDerived>>
    : std::true_type
  {
      /// @brief True for mass forms.
      static constexpr bool Value = true;
  };
}

namespace Rodin::Variational
{
  /**
   * @brief Discrete mass form @f$ m(u,v) = \int_\Omega c\, u \cdot v\,dx @f$.
   *
   * The coefficient @f$ c @f$ is optional: with @p CoefficientDerived
   * @c void the form is the unweighted mass form. It computes the same
   * operator as <tt>Integral(c * u, v)</tt>, respectively
   * <tt>Integral(u, v)</tt>, including the quadrature order that path picks.
   *
   * @tparam CoefficientDerived Derived type of the scalar coefficient, or
   * @c void for none.
   */
  template <class Solution, class TrialFES, class TestFES, class Operator,
    class CoefficientDerived>
  class MassForm final : public BilinearFormBase<Operator>
  {
      using TrialMeshType = typename FormLanguage::Traits<TrialFES>::MeshType;
      using TestMeshType = typename FormLanguage::Traits<TestFES>::MeshType;
      using TrialContextType = typename FormLanguage::Traits<TrialMeshType>::ContextType;
      using TestContextType = typename FormLanguage::Traits<TestMeshType>::ContextType;

    public:
      /// @brief Solution type of the trial function.
      using SolutionType = Solution;

      /// @brief Trial finite element space type.
      using TrialFESType = TrialFES;

      /// @brief Test finite element space type.
      using TestFESType = TestFES;

      /// @brief Assembled operator type.
      using OperatorType = Operator;

      /// @brief Scalar type of the operator entries.
      using ScalarType = typename FormLanguage::Traits<OperatorType>::ScalarType;

      /// @brief Whether the form carries a coefficient.
      static constexpr bool HasCoefficient = !std::is_void_v<CoefficientDerived>;

      /// @brief Coefficient type, meaningful only when @c HasCoefficient.
      using CoefficientType = FunctionBase<CoefficientDerived>;

      /// @brief Owning handle on the coefficient, empty without one.
      using CoefficientPointer = std::conditional_t<HasCoefficient,
        std::unique_ptr<CoefficientType>, std::nullptr_t>;

      /**
       * @brief Local kernel of the mass form.
       *
       * Evaluates the reference bases once per polytope geometry and
       * quadrature order and reuses them across polytopes; only the
       * quadrature weights, the distortion of the map and the coefficient
       * change from one polytope to the next. The kernel owns a copy of the
       * coefficient, so that each thread assembling with its own copy of the
       * kernel also evaluates its own copy of the coefficient.
       */
      class Kernel
      {
        public:
          /// @brief Local matrix type.
          using MatrixType = Math::Matrix<ScalarType>;

          /**
           * @brief Constructs the kernel of the unweighted mass form.
           * @param[in] trialFES Trial space.
           * @param[in] testFES Test space.
           */
          Kernel(const TrialFES& trialFES, const TestFES& testFES)
            requires(!HasCoefficient)
            : m_trialFES(trialFES),
              m_testFES(testFES)
          {}

          /**
           * @brief Constructs the kernel of the weighted mass form.
           * @param[in] trialFES Trial space.
           * @param[in] testFES Test space.
           * @param[in] coefficient Coefficient, copied.
           */
          Kernel(const TrialFES& trialFES, const TestFES& testFES,
            const CoefficientType& coefficient)
            requires HasCoefficient
            : m_trialFES(trialFES),
              m_testFES(testFES),
              m_coefficient(coefficient.copy())
          {}

          /**
           * @brief Copies the kernel, including its own copy of the coefficient.
           * @param[in] other Kernel to copy.
           */
          Kernel(const Kernel& other)
            : m_trialFES(other.m_trialFES),
              m_testFES(other.m_testFES),
              m_coefficient(clone(other.m_coefficient))
          {}

          /**
           * @brief Computes the local mass matrix of @p polytope.
           * @param[out] out Local matrix, test rows by trial columns.
           * @param[in] polytope Polytope to integrate over.
           */
          void compute(MatrixType& out, const Geometry::Polytope& polytope) const
          {
            const size_t d = polytope.getDimension();
            const Index i = polytope.getIndex();
            const auto& trialFE = m_trialFES.get().getFiniteElement(d, i);
            const auto& testFE = m_testFES.get().getFiniteElement(d, i);
            const size_t base = trialFE.getOrder() + testFE.getOrder();
            // The same order Integral(c * u, v) picks: the integrand's order
            // when the coefficient knows its own, the bases' otherwise.
            size_t order = base;
            if constexpr (HasCoefficient)
            {
              if (const auto k = m_coefficient->getOrder(polytope))
                order = *k + base;
            }
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
              ScalarType weight =
                static_cast<ScalarType>(qf.getWeight(qp) * point.getDistortion());
              if constexpr (HasCoefficient)
              {
                weight *= static_cast<ScalarType>(
                  m_coefficient->getValue(IntegrationPoint(point, &qf, qp)));
              }
              out += weight * m_reference[qp];
            }
          }

        private:
          std::reference_wrapper<const TrialFES> m_trialFES;
          std::reference_wrapper<const TestFES> m_testFES;
          CoefficientPointer m_coefficient;
          mutable bool m_cached = false;
          mutable Geometry::Polytope::Type m_geometry = Geometry::Polytope::Type::Point;
          mutable size_t m_order = 0;
          mutable size_t m_trialCount = 0;
          mutable size_t m_testCount = 0;
          mutable std::vector<MatrixType> m_reference;
      };

      /// @brief Local kernel type, as the assembly asks for it.
      using KernelType = Kernel;

      /// @brief Parent class type.
      using Parent = BilinearFormBase<OperatorType>;

      /// @brief Assembly type selected for the trial and test contexts.
      using AssemblyType = typename Assembly::Default<TrialContextType,
        TestContextType>::template Type<OperatorType, MassForm>;

      /**
       * @brief Constructs and assembles the unweighted mass form.
       * @param[in] u Trial function.
       * @param[in] v Test function.
       */
      MassForm(const TrialFunction<SolutionType, TrialFESType>& u,
        const TestFunction<TestFESType>& v)
        requires(!HasCoefficient)
        : m_u(u),
          m_v(v)
      {
        assemble();
      }

      /**
       * @brief Constructs and assembles the weighted mass form.
       * @param[in] c Coefficient, copied.
       * @param[in] u Trial function.
       * @param[in] v Test function.
       */
      MassForm(const CoefficientType& c,
        const TrialFunction<SolutionType, TrialFESType>& u,
        const TestFunction<TestFESType>& v)
        requires HasCoefficient
        : m_u(u),
          m_v(v),
          m_coefficient(c.copy())
      {
        assemble();
      }

      /**
       * @brief Constructs and assembles the mass form weighted by a value
       * lifted to a RealFunction.
       * @param[in] c Constant or callable coefficient.
       * @param[in] u Trial function.
       * @param[in] v Test function.
       */
      template <class L>
        requires std::is_same_v<CoefficientDerived,
          typename FormLanguage::FunctionDerived<RealFunction<L>>::Type>
      MassForm(const L& c, const TrialFunction<SolutionType, TrialFESType>& u,
        const TestFunction<TestFESType>& v)
        : MassForm(RealFunction<L>(c), u, v)
      {}

      /**
       * @brief Copy constructor, copying the coefficient.
       * @param[in] other Form to copy.
       */
      MassForm(const MassForm& other)
        : Parent(other),
          m_u(other.m_u),
          m_v(other.m_v),
          m_coefficient(clone(other.m_coefficient)),
          m_attributes(other.m_attributes),
          m_operator(other.m_operator),
          m_assembly(other.m_assembly)
      {}

      /// @brief Move constructor.
      MassForm(MassForm&&) = default;

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

      /**
       * @brief Gets the coefficient.
       * @pre The form carries a coefficient.
       */
      const CoefficientType& getCoefficient() const
        requires HasCoefficient
      {
        return *m_coefficient;
      }

      /// @brief Gets the region the form integrates over.
      Geometry::Region getRegion() const
      {
        return Geometry::Region::Cells;
      }

      /// @brief Gets the attributes the form is restricted to, empty for all.
      const FlatSet<Geometry::Attribute>& getAttributes() const
      {
        return m_attributes;
      }

      /**
       * @brief Restricts the form to cells carrying one of @p attributes and
       * reassembles it.
       * @param[in] attributes Non-empty set of attributes.
       * @returns Reference to this form.
       */
      MassForm& over(const FlatSet<Geometry::Attribute>& attributes)
      {
        assert(attributes.size() > 0);
        m_attributes = attributes;
        assemble();
        return *this;
      }

      /**
       * @brief Restricts the form to cells carrying the given attributes and
       * reassembles it.
       * @param[in] a1 First attribute.
       * @param[in] as Further attributes.
       * @returns Reference to this form.
       */
      template <class A1, class... As>
      MassForm& over(const A1& a1, const As&... as)
      {
        return over(FlatSet<Geometry::Attribute>{a1, as...});
      }

      /**
       * @brief Builds the kernel the assembly runs on each polytope.
       * @returns Kernel of this form.
       */
      Kernel getKernel() const
      {
        if constexpr (HasCoefficient)
          return Kernel(getTrialFunction().getFiniteElementSpace(),
            getTestFunction().getFiniteElementSpace(), *m_coefficient);
        else
          return Kernel(getTrialFunction().getFiniteElementSpace(),
            getTestFunction().getFiniteElementSpace());
      }

      MassForm* copy() const noexcept override
      {
        return new MassForm(*this);
      }

    private:
      /// @brief Deep-copies a coefficient handle.
      static CoefficientPointer clone(const CoefficientPointer& coefficient)
      {
        if constexpr (HasCoefficient)
          return coefficient ? CoefficientPointer(coefficient->copy()) : nullptr;
        else
          return nullptr;
      }

      std::reference_wrapper<const TrialFunction<SolutionType, TrialFESType>> m_u;
      std::reference_wrapper<const TestFunction<TestFESType>> m_v;
      CoefficientPointer m_coefficient;
      FlatSet<Geometry::Attribute> m_attributes;
      OperatorType m_operator;
      AssemblyType m_assembly;
  };

  /// @brief Deduction guide for the unweighted mass form.
  template <class Solution, class TrialFES, class TestFES>
  MassForm(const TrialFunction<Solution, TrialFES>&, const TestFunction<TestFES>&)
    -> MassForm<Solution, TrialFES, TestFES,
      Math::SparseMatrix<
        typename FormLanguage::Mult<typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>>;

  /// @brief Deduction guide for the mass form weighted by a function.
  template <class CoefficientDerived, class Solution, class TrialFES, class TestFES>
  MassForm(const FunctionBase<CoefficientDerived>&,
    const TrialFunction<Solution, TrialFES>&, const TestFunction<TestFES>&)
    -> MassForm<Solution, TrialFES, TestFES,
      Math::SparseMatrix<
        typename FormLanguage::Mult<typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
      CoefficientDerived>;

  /// @brief Deduction guide for the mass form weighted by a lifted value.
  template <class L, class Solution, class TrialFES, class TestFES>
    requires(!std::is_base_of_v<FormLanguage::Base, L>)
  MassForm(
    const L&, const TrialFunction<Solution, TrialFES>&, const TestFunction<TestFES>&)
    -> MassForm<Solution, TrialFES, TestFES,
      Math::SparseMatrix<
        typename FormLanguage::Mult<typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
      typename FormLanguage::FunctionDerived<RealFunction<L>>::Type>;
}

#endif
