/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BILINEARFORM_H
#define RODIN_VARIATIONAL_BILINEARFORM_H

#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/FormLanguage/List.h"

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
      using NativeAssembly = Assembly::Native<BilinearFormBase>;
      using OpenMPAssembly = Assembly::OpenMP<BilinearFormBase>;

      BilinearFormBase()
      {
        m_assembly.reset(new NativeAssembly);
      }

      BilinearFormBase(const BilinearFormBase& other)
        : FormLanguage::Base(other),
          m_assembly(other.m_assembly->copy()),
          m_bfis(other.m_bfis)
      {}

      BilinearFormBase(BilinearFormBase&& other)
        :  FormLanguage::Base(std::move(other)),
          m_assembly(std::move(other.m_assembly)),
          m_bfis(std::move(other.m_bfis))
      {}

      constexpr
      const FormLanguage::List<BilinearFormIntegratorBase>& getIntegrators() const
      {
        return m_bfis;
      }

      BilinearFormBase& setAssembly(const Assembly::AssemblyBase<BilinearFormBase>& assembly)
      {
        m_assembly.reset(assembly.copy());
        return *this;
      }

      const Assembly::AssemblyBase<BilinearFormBase>& getAssembly() const
      {
        assert(m_assembly);
        return *m_assembly;
      }

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
       * @brief Builds the bilinear form the given bilinear integrator
       * @param[in] bfi Bilinear integrator which will be used to
       * build the bilinear form.
       * @returns Reference to this (for method chaining)
       */
      virtual BilinearFormBase& from(const BilinearFormIntegratorBase& bfi)
      {
        m_bfis.clear();
        add(bfi).assemble();
        return *this;
      }

      virtual BilinearFormBase& from(
          const FormLanguage::List<BilinearFormIntegratorBase>& bfi)
      {
        m_bfis.clear();
        add(bfi).assemble();
        return *this;
      }

      virtual BilinearFormBase& operator=(const BilinearFormIntegratorBase& bfi)
      {
        from(bfi).assemble();
        return *this;
      }

      /**
       * @todo
       */
      virtual BilinearFormBase& operator=(
          const FormLanguage::List<BilinearFormIntegratorBase>& bfis)
      {
        from(bfis).assemble();
        return *this;
      }

      /**
       * @brief Adds a bilinear integrator to the bilinear form.
       * @returns Reference to this (for method chaining)
       */
      virtual BilinearFormBase& add(const BilinearFormIntegratorBase& bfi)
      {
        m_bfis.add(bfi);
        return *this;
      }

      virtual BilinearFormBase& add(const FormLanguage::List<BilinearFormIntegratorBase>& bfis)
      {
        m_bfis.add(bfis);
        return *this;
      }

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

    private:
      std::unique_ptr<Assembly::AssemblyBase<BilinearFormBase>> m_assembly;
      FormLanguage::List<BilinearFormIntegratorBase> m_bfis;
  };

  /**
   * @ingroup BilinearFormSpecializations
   * @brief Speciallization of BilinearForm for Math::SparseMatrix.
   *
   * This specialization aids in the construction of a @f$ n \times m @f$
   * matrix @f$ A @f$, which is associated to the bilinear form. Here, @f$ n
   * @f$ represents the size (total number of degrees-of-freedom) of the trial
   * space, and @f$ m @f$ represents the size of the test space.
   */
  template <class TrialFES, class TestFES>
  class BilinearForm<TrialFES, TestFES, Context::Serial, Math::SparseMatrix> final
    : public BilinearFormBase<Math::SparseMatrix>
  {
    static_assert(
        std::is_same_v<TrialFES, TestFES>,
        "Different trial and test spaces are currently not supported.");

    static_assert(std::is_same_v<typename TrialFES::Context, Context::Serial>);

    public:
      /// Context of BilinearForm
      using Context = typename TrialFES::Context;

      /// Type of operator associated to the bilinear form
      using OperatorType = Math::SparseMatrix;

      /// Parent class
      using Parent = BilinearFormBase<Math::SparseMatrix>;

      /**
       * @brief Constructs a BilinearForm from a TrialFunction and
       * TestFunction.
       *
       * @param[in] u Trial function argument
       * @param[in] v Test function argument
       */
      constexpr
      BilinearForm(const TrialFunction<TrialFES>& u, const TestFunction<TestFES>& v)
        :  m_u(u), m_v(v)
      {}

      constexpr
      BilinearForm(const BilinearForm& other)
        : Parent(other),
          m_u(other.m_u), m_v(other.m_v)
      {}

      constexpr
      BilinearForm(BilinearForm&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)), m_v(std::move(other.m_v)),
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
        assert(false);
        return 0;
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

      BilinearForm& operator=(const BilinearFormIntegratorBase& bfi) override
      {
        from(bfi).assemble();
        return *this;
      }

      /**
       * @todo
       */
      BilinearForm& operator=(
          const FormLanguage::List<BilinearFormIntegratorBase>& bfis) override
      {
        from(bfis).assemble();
        return *this;
      }

      /**
       * @brief Gets the reference to sparse matrix.
       * @returns Reference to the associated sparse matrix.
       */
      virtual OperatorType& getOperator() override
      {
        return m_operator;
      }

      /**
       * @brief Gets the reference to the (local) associated sparse matrix
       * to the bilinear form.
       * @returns Constant reference to the associated sparse matrix.
       */
      virtual const OperatorType& getOperator() const override
      {
        return m_operator;
      }

      virtual BilinearForm* copy() const noexcept override
      {
        return new BilinearForm(*this);
      }

    private:
      std::reference_wrapper<const TrialFunction<TrialFES>> m_u;
      std::reference_wrapper<const TestFunction<TestFES>>   m_v;
      OperatorType m_operator;
  };

  template <class TrialFES, class TestFES>
  BilinearForm(TrialFunction<TrialFES>&, TestFunction<TestFES>&)
    -> BilinearForm<TrialFES, TestFES, typename TrialFES::Context, Math::SparseMatrix>;
}

#include "BilinearForm.hpp"

#endif

