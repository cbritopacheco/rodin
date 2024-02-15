/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BILINEARFORM_H
#define RODIN_VARIATIONAL_BILINEARFORM_H

#include "Rodin/Configure.h"

#include "Rodin/Pair.h"
#include "Rodin/FormLanguage/List.h"
#include "Rodin/Math/SparseMatrix.h"

#include "Rodin/Assembly/ForwardDecls.h"

#include "ForwardDecls.h"
#include "TrialFunction.h"
#include "TestFunction.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
  /**
   * @defgroup BilinearFormSpecializations BilinearForm Template Specializations
   * @brief Template specializations of the BilinearForm class.
   * @see BilinearForm
   */

  template <class OperatorType>
  class BilinearFormBase : public FormLanguage::Base
  {
    public:
      BilinearFormBase()
      {}

      BilinearFormBase(const BilinearFormBase& other)
        : FormLanguage::Base(other)
      {}

      BilinearFormBase(BilinearFormBase&& other)
        : FormLanguage::Base(std::move(other))
      {}

      /**
       * @brief Assembles the bilinear form.
       *
       * This method will assemble the underlying sparse matrix associated
       * the bilinear form.
       *
       * @see getMatrix()
       */
      virtual void assemble() = 0;

      /**
       * @brief Gets the reference to the associated operator of the bilinear
       * form.
       */
      virtual OperatorType& getOperator() = 0;

      /** @brief Gets a constant reference to the associated operator of the
       * bilinear form.
       */
      virtual const OperatorType& getOperator() const = 0;

      /**
       * @brief Gets the reference to the associated TrialFunction object.
       * @returns Reference to this (for method chaining)
       */
      virtual const FormLanguage::Base& getTrialFunction() const = 0;

      /**
       * @brief Gets the reference to the associated TestFunction object.
       * @returns Reference to this (for method chaining)
       */
      virtual const FormLanguage::Base& getTestFunction() const = 0;

      virtual BilinearFormBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup BilinearFormSpecializations
   * @brief Speciallization of BilinearForm for a matrix type.
   *
   * This specialization aids in the construction of a @f$ m \times n @f$
   * matrix @f$ A @f$, which is associated to the bilinear form. Here, @f$ n
   * @f$ represents the size (total number of degrees-of-freedom) of the trial
   * space, and @f$ m @f$ represents the size of the test space.
   */
  template <class TrialFES, class TestFES, class MatrixType>
  class BilinearForm final
    : public BilinearFormBase<MatrixType>
  {
    public:
      /// Type of operator associated to the bilinear form
      using OperatorType = MatrixType;

      /// Parent class
      using Parent = BilinearFormBase<MatrixType>;

      using SequentialAssembly = Assembly::Sequential<OperatorType, BilinearForm>;

      using MultithreadedAssembly = Assembly::Multithreaded<OperatorType, BilinearForm>;

      constexpr
      BilinearForm(Pair<
          std::reference_wrapper<const TrialFunction<TrialFES>>,
          std::reference_wrapper<const TestFunction<TestFES>>> uv)
        : BilinearForm(uv.first(), uv.second())
      {}

      /**
       * @brief Constructs a BilinearForm from a TrialFunction and
       * TestFunction.
       *
       * @param[in] u Trial function argument
       * @param[in] v Test function argument
       */
      constexpr
      BilinearForm(
          std::reference_wrapper<const TrialFunction<TrialFES>> u,
          std::reference_wrapper<const TestFunction<TestFES>> v)
        :  m_u(u), m_v(v)
      {
#ifdef RODIN_MULTITHREADED
        m_assembly.reset(new MultithreadedAssembly);
#else
        m_assembly.reset(new SequentialAssembly);
#endif
      }

      constexpr
      BilinearForm(const BilinearForm& other)
        : Parent(other),
          m_u(other.m_u), m_v(other.m_v),
          m_assembly(other.m_assembly->copy()),
          m_lbfis(other.m_lbfis),
          m_gbfis(other.m_gbfis)
      {}

      constexpr
      BilinearForm(BilinearForm&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)), m_v(std::move(other.m_v)),
          m_assembly(std::move(other.m_assembly)),
          m_lbfis(std::move(other.m_lbfis)),
          m_gbfis(std::move(other.m_gbfis)),
          m_operator(std::move(other.m_operator))
      {}

      /**
       * @brief Evaluates the linear form at the functions @f$ u @f$ and @f$
       * v @f$.
       *
       * Given grid functions @f$ u @f$ and @f$ v @f$, this function will
       * compute the action of the bilinear mapping @f$ a(u, v) @f$.
       *
       * @returns The action @f$ a(u, v) @f$ which the bilinear form takes
       * at @f$ ( u, v ) @f$.
       */
      constexpr
      Scalar operator()(const GridFunction<TrialFES>& u, const GridFunction<TestFES>& v) const
      {
        const auto& trialWeights = u.getWeights();
        const auto& testWeights = v.getWeights();
        if (!trialWeights.has_value())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Trial GridFunction weights have not been calculated. "
            << "Call " << Alert::Identifier::Function("setWeights()")
            << " on the GridFunction object."
            << Alert::Raise;
        }
        assert(trialWeights.has_value());
        if (!testWeights.has_value())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Test GridFunction weights have not been calculated. "
            << "Call " << Alert::Identifier::Function("setWeights()")
            << " on the GridFunction object."
            << Alert::Raise;
        }
        assert(testWeights.has_value());
        return (getOperator() * testWeights.value()).dot(trialWeights.value());
      }

      void assemble() override;

      const TrialFunction<TrialFES>& getTrialFunction() const override
      {
        return m_u.get();
      }

      const TestFunction<TestFES>& getTestFunction() const override
      {
        return m_v.get();
      }

      BilinearForm& operator=(const LocalBilinearFormIntegratorBase& bfi)
      {
        this->from(bfi).assemble();
        return *this;
      }

      /**
       * @todo
       */
      BilinearForm& operator=(
          const FormLanguage::List<LocalBilinearFormIntegratorBase>& bfis)
      {
        this->from(bfis).assemble();
        return *this;
      }

      BilinearForm& operator=(const GlobalBilinearFormIntegratorBase& bfi)
      {
        this->from(bfi).assemble();
        return *this;
      }

      /**
       * @todo
       */
      BilinearForm& operator=(
          const FormLanguage::List<GlobalBilinearFormIntegratorBase>& bfis)
      {
        this->from(bfis).assemble();
        return *this;
      }

      constexpr
      FormLanguage::List<LocalBilinearFormIntegratorBase>& getLocalIntegrators()
      {
        return m_lbfis;
      }

      constexpr
      const FormLanguage::List<LocalBilinearFormIntegratorBase>& getLocalIntegrators() const
      {
        return m_lbfis;
      }

      constexpr
      FormLanguage::List<GlobalBilinearFormIntegratorBase>& getGlobalIntegrators()
      {
        return m_gbfis;
      }

      constexpr
      const FormLanguage::List<GlobalBilinearFormIntegratorBase>& getGlobalIntegrators() const
      {
        return m_gbfis;
      }

      /**
       * @brief Gets the reference to sparse matrix.
       * @returns Reference to the associated sparse matrix.
       */
      OperatorType& getOperator() override
      {
        return m_operator;
      }

      /**
       * @brief Gets the reference to the (local) associated sparse matrix
       * to the bilinear form.
       * @returns Constant reference to the associated sparse matrix.
       */
      const OperatorType& getOperator() const override
      {
        return m_operator;
      }

      BilinearForm& setAssembly(const Assembly::AssemblyBase<OperatorType, BilinearForm>& assembly)
      {
        m_assembly.reset(assembly.copy());
        return *this;
      }

      const Assembly::AssemblyBase<OperatorType, BilinearForm>& getAssembly() const
      {
        assert(m_assembly);
        return *m_assembly;
      }

      /**
       * @brief Builds the bilinear form the given bilinear integrator
       * @param[in] bfi Bilinear integrator which will be used to
       * build the bilinear form.
       * @returns Reference to this (for method chaining)
       */
      virtual BilinearForm& from(const LocalBilinearFormIntegratorBase& bfi)
      {
        m_lbfis.clear();
        add(bfi).assemble();
        return *this;
      }

      virtual BilinearForm& from(
          const FormLanguage::List<LocalBilinearFormIntegratorBase>& bfi)
      {
        m_lbfis.clear();
        add(bfi).assemble();
        return *this;
      }

      /**
       * @brief Adds a bilinear integrator to the bilinear form.
       * @returns Reference to this (for method chaining)
       */
      virtual BilinearForm& add(const LocalBilinearFormIntegratorBase& bfi)
      {
        m_lbfis.add(bfi);
        return *this;
      }

      virtual BilinearForm& add(const FormLanguage::List<LocalBilinearFormIntegratorBase>& bfis)
      {
        m_lbfis.add(bfis);
        return *this;
      }

      /**
       * @brief Builds the bilinear form the given bilinear integrator
       * @param[in] bfi Bilinear integrator which will be used to
       * build the bilinear form.
       * @returns Reference to this (for method chaining)
       */
      virtual BilinearForm& from(const GlobalBilinearFormIntegratorBase& bfi)
      {
        m_gbfis.clear();
        add(bfi).assemble();
        return *this;
      }

      virtual BilinearForm& from(
          const FormLanguage::List<GlobalBilinearFormIntegratorBase>& bfi)
      {
        m_gbfis.clear();
        add(bfi).assemble();
        return *this;
      }

      /**
       * @brief Adds a bilinear integrator to the bilinear form.
       * @returns Reference to this (for method chaining)
       */
      virtual BilinearForm& add(const GlobalBilinearFormIntegratorBase& bfi)
      {
        m_gbfis.add(bfi);
        return *this;
      }

      virtual BilinearForm& add(const FormLanguage::List<GlobalBilinearFormIntegratorBase>& bfis)
      {
        m_gbfis.add(bfis);
        return *this;
      }

      inline BilinearForm* copy() const noexcept override
      {
        return new BilinearForm(*this);
      }

    private:
      std::reference_wrapper<const TrialFunction<TrialFES>> m_u;
      std::reference_wrapper<const TestFunction<TestFES>>   m_v;
      std::unique_ptr<Assembly::AssemblyBase<OperatorType, BilinearForm>> m_assembly;
      FormLanguage::List<LocalBilinearFormIntegratorBase>        m_lbfis;
      FormLanguage::List<GlobalBilinearFormIntegratorBase>       m_gbfis;
      OperatorType m_operator;
  };

  template <class TrialFES, class TestFES>
  BilinearForm(TrialFunction<TrialFES>&, TestFunction<TestFES>&)
    -> BilinearForm<TrialFES, TestFES, Math::SparseMatrix>;
}

#include "BilinearForm.hpp"

#endif

