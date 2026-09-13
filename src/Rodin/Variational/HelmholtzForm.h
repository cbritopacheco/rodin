/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_HELMHOLTZFORM_H
#define RODIN_VARIATIONAL_HELMHOLTZFORM_H

#include <functional>
#include <memory>
#include <type_traits>

#include "Rodin/Assembly/Default.h"
#include "Rodin/FormLanguage/Traits.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Geometry/Region.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SparseMatrix.h"

#include "BilinearForm.h"
#include "DiffusionForm.h"
#include "Function.h"
#include "MassForm.h"
#include "RealFunction.h"

namespace Rodin::FormLanguage
{
  /// @brief Traits of a Variational::HelmholtzForm.
  template <class Solution, class TrialFES, class TestFES, class Operator,
    class DiffusionDerived, class MassDerived>
  struct Traits<Variational::HelmholtzForm<Solution, TrialFES, TestFES, Operator,
    DiffusionDerived, MassDerived>>
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

  /// @brief Marks Variational::HelmholtzForm as a named form.
  template <class Solution, class TrialFES, class TestFES, class Operator,
    class DiffusionDerived, class MassDerived>
  struct IsNamedForm<Variational::HelmholtzForm<Solution, TrialFES, TestFES, Operator,
    DiffusionDerived, MassDerived>> : std::true_type
  {
      /// @brief True for Helmholtz forms.
      static constexpr bool Value = true;
  };
}

namespace Rodin::Variational
{
  /**
   * @brief Discrete Helmholtz-type form
   * @f$ h(u,v) = \int_\Omega a\, \nabla u \cdot \nabla v\,dx
   *   + \int_\Omega c\, u\, v\,dx @f$.
   *
   * Computes the same operator as
   * <tt>Integral(a * Grad(u), Grad(v)) + Integral(c * u, v)</tt>. Each term is
   * integrated at the quadrature order its own integral would pick, so the
   * two terms may use different quadratures on the same polytope; what the
   * form saves over the sum of two forms is the second pass over the mesh and
   * the second scatter, not a quadrature. The sign of @f$ c @f$ is the
   * caller's: a negative @f$ c = -k^2 @f$ gives the indefinite Helmholtz
   * operator, a positive one the reaction-diffusion operator used by the
   * Hilbertian regularizations. Scalar spaces only.
   *
   * @tparam DiffusionDerived Derived type of the diffusion coefficient.
   * @tparam MassDerived Derived type of the mass coefficient.
   */
  template <class Solution, class TrialFES, class TestFES, class Operator,
    class DiffusionDerived, class MassDerived>
  class HelmholtzForm final : public BilinearFormBase<Operator>
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

      /// @brief Diffusion coefficient type.
      using DiffusionCoefficientType = FunctionBase<DiffusionDerived>;

      /// @brief Mass coefficient type.
      using MassCoefficientType = FunctionBase<MassDerived>;

      /// @brief Kernel of the diffusion term.
      using DiffusionKernelType = typename DiffusionForm<Solution, TrialFES, TestFES,
        Operator, DiffusionDerived>::KernelType;

      /// @brief Kernel of the mass term.
      using MassKernelType =
        typename MassForm<Solution, TrialFES, TestFES, Operator, MassDerived>::KernelType;

      /**
       * @brief Local kernel of the Helmholtz form.
       *
       * Runs the diffusion and the mass kernels on the same polytope and sums
       * their local matrices.
       */
      class Kernel
      {
        public:
          /// @brief Local matrix type.
          using MatrixType = Math::Matrix<ScalarType>;

          /**
           * @brief Constructs the kernel.
           * @param[in] trialFES Trial space.
           * @param[in] testFES Test space.
           * @param[in] a Diffusion coefficient, copied.
           * @param[in] c Mass coefficient, copied.
           */
          Kernel(const TrialFES& trialFES, const TestFES& testFES,
            const DiffusionCoefficientType& a, const MassCoefficientType& c)
            : m_diffusion(trialFES, testFES, a),
              m_mass(trialFES, testFES, c)
          {}

          /**
           * @brief Computes the local Helmholtz matrix of @p polytope.
           * @param[out] out Local matrix, test rows by trial columns.
           * @param[in] polytope Polytope to integrate over.
           */
          void compute(MatrixType& out, const Geometry::Polytope& polytope) const
          {
            m_diffusion.compute(out, polytope);
            m_mass.compute(m_scratch, polytope);
            out += m_scratch;
          }

        private:
          DiffusionKernelType m_diffusion;
          MassKernelType m_mass;
          mutable MatrixType m_scratch;
      };

      /// @brief Local kernel type, as the assembly asks for it.
      using KernelType = Kernel;

      /// @brief Parent class type.
      using Parent = BilinearFormBase<OperatorType>;

      /// @brief Assembly type selected for the trial and test contexts.
      using AssemblyType = typename Assembly::Default<TrialContextType,
        TestContextType>::template Type<OperatorType, HelmholtzForm>;

      /**
       * @brief Constructs and assembles the Helmholtz form.
       *
       * Each coefficient is either a function, copied, or a value lifted to a
       * RealFunction.
       *
       * @param[in] a Diffusion coefficient.
       * @param[in] c Mass coefficient.
       * @param[in] u Trial function.
       * @param[in] v Test function.
       */
      template <class A, class C>
      HelmholtzForm(const A& a, const C& c,
        const TrialFunction<SolutionType, TrialFESType>& u,
        const TestFunction<TestFESType>& v)
        : m_u(u),
          m_v(v),
          m_diffusionCoefficient(lift<DiffusionDerived>(a)),
          m_massCoefficient(lift<MassDerived>(c))
      {
        assemble();
      }

      /**
       * @brief Copy constructor, copying the coefficients.
       * @param[in] other Form to copy.
       */
      HelmholtzForm(const HelmholtzForm& other)
        : Parent(other),
          m_u(other.m_u),
          m_v(other.m_v),
          m_diffusionCoefficient(other.m_diffusionCoefficient->copy()),
          m_massCoefficient(other.m_massCoefficient->copy()),
          m_attributes(other.m_attributes),
          m_operator(other.m_operator),
          m_assembly(other.m_assembly)
      {}

      /// @brief Move constructor.
      HelmholtzForm(HelmholtzForm&&) = default;

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

      /// @brief Gets the diffusion coefficient.
      const DiffusionCoefficientType& getDiffusionCoefficient() const
      {
        return *m_diffusionCoefficient;
      }

      /// @brief Gets the mass coefficient.
      const MassCoefficientType& getMassCoefficient() const
      {
        return *m_massCoefficient;
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
      HelmholtzForm& over(const FlatSet<Geometry::Attribute>& attributes)
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
      HelmholtzForm& over(const A1& a1, const As&... as)
      {
        return over(FlatSet<Geometry::Attribute>{a1, as...});
      }

      /**
       * @brief Builds the kernel the assembly runs on each polytope.
       * @returns Kernel of this form.
       */
      Kernel getKernel() const
      {
        return Kernel(getTrialFunction().getFiniteElementSpace(),
          getTestFunction().getFiniteElementSpace(), *m_diffusionCoefficient,
          *m_massCoefficient);
      }

      HelmholtzForm* copy() const noexcept override
      {
        return new HelmholtzForm(*this);
      }

    private:
      /**
       * @brief Copies @p value as a coefficient, lifting a non-function to a
       * RealFunction.
       */
      template <class Derived, class Value>
      static std::unique_ptr<FunctionBase<Derived>> lift(const Value& value)
      {
        if constexpr (std::is_base_of_v<FunctionBase<Derived>, Value>)
        {
          return std::unique_ptr<FunctionBase<Derived>>(
            static_cast<const FunctionBase<Derived>&>(value).copy());
        }
        else
        {
          static_assert(std::is_same_v<Derived,
            typename FormLanguage::FunctionDerived<RealFunction<Value>>::Type>);
          return std::unique_ptr<FunctionBase<Derived>>(
            static_cast<FunctionBase<Derived>*>(RealFunction<Value>(value).copy()));
        }
      }

      std::reference_wrapper<const TrialFunction<SolutionType, TrialFESType>> m_u;
      std::reference_wrapper<const TestFunction<TestFESType>> m_v;
      std::unique_ptr<DiffusionCoefficientType> m_diffusionCoefficient;
      std::unique_ptr<MassCoefficientType> m_massCoefficient;
      FlatSet<Geometry::Attribute> m_attributes;
      OperatorType m_operator;
      AssemblyType m_assembly;
  };

  /// @brief Deduction guide for two function coefficients.
  template <class DiffusionDerived, class MassDerived, class Solution, class TrialFES,
    class TestFES>
  HelmholtzForm(const FunctionBase<DiffusionDerived>&, const FunctionBase<MassDerived>&,
    const TrialFunction<Solution, TrialFES>&, const TestFunction<TestFES>&)
    -> HelmholtzForm<Solution, TrialFES, TestFES,
      Math::SparseMatrix<
        typename FormLanguage::Mult<typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
      DiffusionDerived, MassDerived>;

  /// @brief Deduction guide for a function diffusion and a lifted mass coefficient.
  template <class DiffusionDerived, class C, class Solution, class TrialFES,
    class TestFES>
    requires(!std::is_base_of_v<FormLanguage::Base, C>)
  HelmholtzForm(const FunctionBase<DiffusionDerived>&, const C&,
    const TrialFunction<Solution, TrialFES>&, const TestFunction<TestFES>&)
    -> HelmholtzForm<Solution, TrialFES, TestFES,
      Math::SparseMatrix<
        typename FormLanguage::Mult<typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
      DiffusionDerived, typename FormLanguage::FunctionDerived<RealFunction<C>>::Type>;

  /// @brief Deduction guide for a lifted diffusion and a function mass coefficient.
  template <class A, class MassDerived, class Solution, class TrialFES, class TestFES>
    requires(!std::is_base_of_v<FormLanguage::Base, A>)
  HelmholtzForm(const A&, const FunctionBase<MassDerived>&,
    const TrialFunction<Solution, TrialFES>&, const TestFunction<TestFES>&)
    -> HelmholtzForm<Solution, TrialFES, TestFES,
      Math::SparseMatrix<
        typename FormLanguage::Mult<typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
      typename FormLanguage::FunctionDerived<RealFunction<A>>::Type, MassDerived>;

  /// @brief Deduction guide for two lifted coefficients.
  template <class A, class C, class Solution, class TrialFES, class TestFES>
    requires(!std::is_base_of_v<FormLanguage::Base, A> &&
              !std::is_base_of_v<FormLanguage::Base, C>)
  HelmholtzForm(const A&, const C&, const TrialFunction<Solution, TrialFES>&,
    const TestFunction<TestFES>&)
    -> HelmholtzForm<Solution, TrialFES, TestFES,
      Math::SparseMatrix<
        typename FormLanguage::Mult<typename FormLanguage::Traits<TrialFES>::ScalarType,
          typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
      typename FormLanguage::FunctionDerived<RealFunction<A>>::Type,
      typename FormLanguage::FunctionDerived<RealFunction<C>>::Type>;
}

#endif
