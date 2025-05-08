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

#include "Exceptions/TrialFunctionMismatchException.h"
#include "Exceptions/TestFunctionMismatchException.h"

#include "ForwardDecls.h"
#include "TrialFunction.h"
#include "TestFunction.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::FormLanguage
{
  template <class Operator>
  struct Traits<Variational::BilinearFormBase<Operator>>
  {
    using ScalarType = typename FormLanguage::Traits<Operator>::ScalarType;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup BilinearFormSpecializations BilinearForm Template Specializations
   * @brief Template specializations of the BilinearForm class.
   * @see BilinearForm
   */

  template <class Operator>
  class BilinearFormBase : public FormLanguage::Base
  {
    public:
      using OperatorType = Operator;

      using Parent = FormLanguage::Base;

      /**
       * @brief Constructs a linear form with a default constructed vector
       * which is owned by the LinearFormBase instance.
       */
      BilinearFormBase()
        : m_operator(OperatorType())
      {}

      /**
       * @brief Constructs a linear form with reference to vector which is not
       * owned by the LinearFormBase instance.
       */
      BilinearFormBase(OperatorType& vec)
        : m_operator(std::ref(vec))
      {}

      /**
       * @brief Constructs a linear form with a vector which is owned by the
       * LinearFormBase instance.
       */
      BilinearFormBase(OperatorType&& vec)
        : m_operator(std::move(vec))
      {}

      BilinearFormBase(const BilinearFormBase& other)
        : Parent(other),
          m_operator(other.m_operator)
      {}

      BilinearFormBase(BilinearFormBase&& other)
        : Parent(std::move(other)),
          m_operator(std::move(other.m_operator))
      {}

      /**
       * @brief Gets the reference to the associated operator of the bilinear
       * form.
       */
      OperatorType& getOperator()
      {
        auto& ref = std::visit([](auto& m) -> OperatorType& { return m; }, m_operator);
        return ref;
      }

      /** @brief Gets a constant reference to the associated operator of the
       * bilinear form.
       */
      const OperatorType& getOperator() const
      {
        const auto& ref = std::visit([](const auto& m) -> const OperatorType& { return m; }, m_operator);
        return ref;
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
      std::variant<std::reference_wrapper<OperatorType>, OperatorType> m_operator;
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
  template <class TrialFES, class TestFES, class Operator>
  class BilinearForm final
    : public BilinearFormBase<Operator>
  {
    public:
      using TrialFESScalarType  = typename FormLanguage::Traits<TrialFES>::ScalarType;

      using TestFESScalarType   = typename FormLanguage::Traits<TestFES>::ScalarType;

      using ScalarType = typename FormLanguage::Mult<TrialFESScalarType, TestFESScalarType>::Type;

      /// Type of operator associated to the bilinear form
      using OperatorType = Operator;

      using LocalBilinearFormIntegratorBaseType =
        LocalBilinearFormIntegratorBase<ScalarType>;

      using GlobalBilinearFormIntegratorBaseType =
        GlobalBilinearFormIntegratorBase<ScalarType>;

      using LocalBilinearFormIntegratorBaseListType =
        FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      using GlobalBilinearFormIntegratorBaseListType =
        FormLanguage::List<GlobalBilinearFormIntegratorBaseType>;

      /// Parent class
      using Parent = BilinearFormBase<OperatorType>;

      using SequentialAssembly = Assembly::Sequential<OperatorType, BilinearForm>;

      using MultithreadedAssembly = Assembly::Multithreaded<OperatorType, BilinearForm>;

      /**
       * @brief Constructs a LinearForm with a reference to a TestFunction and
       * a default constructed vector owned by the LinearForm instance.
       * @param[in] v Reference to a TestFunction
       */
      constexpr
      BilinearForm(const TrialFunction<TrialFES>& u, const TestFunction<TestFES>& v)
        : BilinearForm(u, v, OperatorType())
      {}

      /**
       * @brief Constructs a LinearForm with a reference to a TestFunction and
       * an non-owned vector.
       * @param[in] v Reference to a TestFunction
       * @param[in] vec Reference to a vector
       */
      constexpr
      BilinearForm(const TrialFunction<TrialFES>& u, const TestFunction<TestFES>& v, OperatorType& op)
        : Parent(op),
          m_u(u), m_v(v)
      {
#ifdef RODIN_MULTITHREADED
        m_assembly.reset(new MultithreadedAssembly);
#else
        m_assembly.reset(new SequentialAssembly);
#endif
      }

      /**
       * @brief Constructs a LinearForm with a references to a TrialFunction and
       * a TestFunction, and an owned operator.
       */
      constexpr
      BilinearForm(const TrialFunction<TrialFES>& u, const TestFunction<TestFES>& v, Operator&& op)
        : Parent(std::move(op)),
          m_u(u), m_v(v)
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
          m_gbfis(std::move(other.m_gbfis))
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
      ScalarType operator()(const GridFunction<TrialFES>& u, const GridFunction<TestFES>& v) const
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
        return (this->getOperator() * testWeights.value()).dot(trialWeights.value());
      }

      void assemble() override
      {
         const auto& trialFES = getTrialFunction().getFiniteElementSpace();
         const auto& testFES = getTestFunction().getFiniteElementSpace();
         getAssembly().execute(this->getOperator(), {
             trialFES, testFES, getLocalIntegrators(), getGlobalIntegrators() });
      }

      const TrialFunction<TrialFES>& getTrialFunction() const override
      {
        return m_u.get();
      }

      const TestFunction<TestFES>& getTestFunction() const override
      {
        return m_v.get();
      }

      BilinearForm& operator=(const LocalBilinearFormIntegratorBaseType& bfi)
      {
        this->from(bfi);
        return *this;
      }

      /**
       * @todo
       */
      BilinearForm& operator=(
          const FormLanguage::List<LocalBilinearFormIntegratorBaseType>& bfis)
      {
        this->from(bfis);
        return *this;
      }

      BilinearForm& operator=(const GlobalBilinearFormIntegratorBaseType& bfi)
      {
        this->from(bfi);
        return *this;
      }

      /**
       * @todo
       */
      BilinearForm& operator=(
          const FormLanguage::List<GlobalBilinearFormIntegratorBaseType>& bfis)
      {
        this->from(bfis);
        return *this;
      }

      constexpr
      LocalBilinearFormIntegratorBaseListType& getLocalIntegrators()
      {
        return m_lbfis;
      }

      constexpr
      const LocalBilinearFormIntegratorBaseListType& getLocalIntegrators() const
      {
        return m_lbfis;
      }

      constexpr
      GlobalBilinearFormIntegratorBaseListType& getGlobalIntegrators()
      {
        return m_gbfis;
      }

      constexpr
      const GlobalBilinearFormIntegratorBaseListType& getGlobalIntegrators() const
      {
        return m_gbfis;
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
      BilinearForm& from(const LocalBilinearFormIntegratorBaseType& bfi)
      {
        m_lbfis.clear();
        add(bfi);
        return *this;
      }

      BilinearForm& from(const LocalBilinearFormIntegratorBaseListType& bfi)
      {
        m_lbfis.clear();
        add(bfi);
        return *this;
      }

      /**
       * @brief Adds a bilinear integrator to the bilinear form.
       * @returns Reference to this (for method chaining)
       */
      BilinearForm& add(const LocalBilinearFormIntegratorBaseType& bfi)
      {
        if (bfi.getTrialFunction().getUUID() != getTrialFunction().getUUID())
          TrialFunctionMismatchException(bfi.getTrialFunction()) << Alert::Raise;
        if (bfi.getTestFunction().getUUID() != getTestFunction().getUUID())
          TestFunctionMismatchException(bfi.getTestFunction()) << Alert::Raise;
        m_lbfis.add(bfi);
        return *this;
      }

      BilinearForm& add(const LocalBilinearFormIntegratorBaseListType& bfis)
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
      BilinearForm& from(const GlobalBilinearFormIntegratorBaseType& bfi)
      {
        m_gbfis.clear();
        add(bfi);
        return *this;
      }

      BilinearForm& from(const GlobalBilinearFormIntegratorBaseListType& bfi)
      {
        m_gbfis.clear();
        add(bfi);
        return *this;
      }

      /**
       * @brief Adds a bilinear integrator to the bilinear form.
       * @returns Reference to this (for method chaining)
       */
      BilinearForm& add(const GlobalBilinearFormIntegratorBaseType& bfi)
      {
        if (bfi.getTrialFunction().getUUID() != getTrialFunction().getUUID())
          TrialFunctionMismatchException(bfi.getTrialFunction()) << Alert::Raise;
        if (bfi.getTestFunction().getUUID() != getTestFunction().getUUID())
          TestFunctionMismatchException(bfi.getTestFunction()) << Alert::Raise;
        m_gbfis.add(bfi);
        return *this;
      }

      BilinearForm& add(const GlobalBilinearFormIntegratorBaseListType& bfis)
      {
        m_gbfis.add(bfis);
        return *this;
      }

      BilinearForm& clear()
      {
        m_lbfis.clear();
        m_gbfis.clear();
        return *this;
      }

      BilinearForm* copy() const noexcept override
      {
        return new BilinearForm(*this);
      }

    private:
      std::reference_wrapper<const TrialFunction<TrialFES>> m_u;
      std::reference_wrapper<const TestFunction<TestFES>>   m_v;
      std::unique_ptr<Assembly::AssemblyBase<OperatorType, BilinearForm>> m_assembly;
      LocalBilinearFormIntegratorBaseListType               m_lbfis;
      GlobalBilinearFormIntegratorBaseListType              m_gbfis;
  };

  template <class TrialFES, class TestFES>
  BilinearForm(const TrialFunction<TrialFES>& u, const TestFunction<TestFES>& v)
    -> BilinearForm<TrialFES, TestFES,
        Math::SparseMatrix<
          typename FormLanguage::Mult<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>
          ::Type>>;

  template <class TrialFES, class TestFES, class Operator>
  BilinearForm(const TrialFunction<TrialFES>& u, const TestFunction<TestFES>& v, Operator&& op)
    -> BilinearForm<TrialFES, TestFES, Operator>;

  template <class TrialFES, class TestFES, class Operator>
  BilinearForm(const TrialFunction<TrialFES>& u, const TestFunction<TestFES>& v, Operator& op)
    -> BilinearForm<TrialFES, TestFES, Operator>;
}

#include "BilinearForm.hpp"

#endif

